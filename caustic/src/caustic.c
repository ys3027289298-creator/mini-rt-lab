
/* caustic: specular-diffuse path-class sampling for caustics in flatland.
   Paths are enumerated by their mirror chain (S events) ending on the
   diffuse floor (D event); class k=0 is direct light (length-1 path).
   Each class is estimated by two strategies (light area sampling and
   receiver cosine sampling through the unfolded specular chain) and
   combined with the power heuristic. A partition-of-unity gate must
   pass before any image is composited. Estimates are checked against a
   high-resolution quadrature reference with Hoeffding/Chebyshev error
   upper bounds. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define EPS 1e-9
#define MAX_MIR 4
#define MAX_K 3
#define MAX_CHAIN 32
#define PIXELS 121
#define XMIN (-3.0)
#define XMAX 3.0

typedef struct { double x, y; } V2;
typedef struct { V2 a, b; double rho; } Mirror;

static V2 vadd(V2 a, V2 b) { V2 r = { a.x + b.x, a.y + b.y }; return r; }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vscale(V2 a, double s) { V2 r = { a.x * s, a.y * s }; return r; }
static double vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static double vcross(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }
static double vlen(V2 a) { return sqrt(vdot(a, a)); }
static V2 vnorm(V2 a) { double l = vlen(a); return l > 0 ? vscale(a, 1.0 / l) : a; }

/* reflect direction d across the infinite line through a,b */
static V2 reflect_dir(V2 d, V2 a, V2 b) {
    V2 nline = vnorm(vsub(b, a));
    /* unit normal of the line */
    V2 n = { -nline.y, nline.x };
    return vsub(d, vscale(n, 2.0 * vdot(d, n)));
}
/* reflect point p across the infinite line through a,b */
static V2 reflect_point(V2 p, V2 a, V2 b) {
    V2 e = vsub(b, a);
    double t = vdot(vsub(p, a), e) / vdot(e, e);
    V2 foot = vadd(a, vscale(e, t));
    return vsub(vscale(foot, 2.0), p);
}

/* ray o+t d vs segment a+u(b-a); t>EPS, u in [0,1] */
static int ray_seg(V2 o, V2 d, V2 a, V2 b, double* t, double* u) {
    V2 e = vsub(b, a);
    double den = vcross(d, e);
    if (fabs(den) < 1e-14) return 0;
    double tt = vcross(vsub(a, o), e) / den;
    double uu = vcross(vsub(a, o), d) / den;
    if (tt <= EPS || uu < -1e-10 || uu > 1.0 + 1e-10) return 0;
    if (t) *t = tt;
    if (u) *u = uu;
    return 1;
}

/* distance from point p to segment a-b */
static double point_seg_dist(V2 p, V2 a, V2 b) {
    V2 e = vsub(b, a);
    double t = vdot(vsub(p, a), e) / vdot(e, e);
    t = t < 0 ? 0 : (t > 1 ? 1 : t);
    return vlen(vsub(p, vadd(a, vscale(e, t))));
}
static double seg_seg_dist(V2 a, V2 b, V2 c, V2 d) {
    double v0 = point_seg_dist(a, c, d), v1 = point_seg_dist(b, c, d);
    double v2 = point_seg_dist(c, a, b), v3 = point_seg_dist(d, a, b);
    return fmin(fmin(v0, v1), fmin(v2, v3));
}

/* ---- scene ---------------------------------------------------------------- */

static V2 g_l0 = { -0.8, 3.0 }, g_l1 = { 0.8, 3.0 };
static double g_Le = 1.0;
static Mirror g_mir[MAX_MIR];
static int g_nmir;
static double g_L; /* light length */

static void scene_build(void) {
    g_L = vlen(vsub(g_l1, g_l0));
    g_nmir = 0;
    /* two near-vertical, slightly tilted mirrors facing each other:
       light can bounce between them, giving real k=2 caustic paths. */
    g_mir[g_nmir++] = (Mirror){ { -1.05, 0.25 }, { -0.95, 2.35 }, 0.9 };
    g_mir[g_nmir++] = (Mirror){ { 1.05, 0.25 }, { 0.95, 2.35 }, 0.9 };
}

/* nearest mirror hit along a ray; returns mirror index or -1 */
static int trace_mirror(V2 o, V2 d, double tmax, double* tout) {
    int best = -1; double bt = 1e30;
    for (int i = 0; i < g_nmir; i++) {
        double t, u;
        if (ray_seg(o, d, g_mir[i].a, g_mir[i].b, &t, &u) && t < bt &&
            (tmax <= 0 || t < tmax)) { bt = t; best = i; }
    }
    if (tout) *tout = bt;
    return best;
}

/* ray vs horizontal light segment; 1 if hit before tmax */
static int trace_light(V2 o, V2 d, double tmax, double* xl, double* tl) {
    if (d.y <= 1e-12) return 0;
    double t = (g_l0.y - o.y) / d.y;
    if (t <= EPS || (tmax > 0 && t >= tmax)) return 0;
    double x = o.x + t * d.x;
    if (x < g_l0.x - 1e-10 || x > g_l1.x + 1e-10) return 0;
    if (xl) *xl = x;
    if (tl) *tl = t;
    return 1;
}

/* ---- path class solver ---------------------------------------------------- */

typedef struct {
    int ok;
    int chain[MAX_K], k;
    double xl;
    double r0, rk, cL, cosX, dtheta;
    double G0, Gk, rho;
} Path;

/* Solve path class `chain` (k specular events) for floor point x and light
   coordinate xl. Uses mirror unfolding and verifies by forward tracing. */
static int solve_chain(double x, const int* chain, int k, double xl, Path* p) {
    memset(p, 0, sizeof(*p));
    p->k = k;
    memcpy(p->chain, chain, sizeof(int) * k);
    p->xl = xl;
    p->rho = 1.0;
    for (int j = 0; j < k; j++) p->rho *= g_mir[chain[j]].rho;

    /* Virtual-image unfolding: reflect the receiver x across the chain to
       get x*. The light stays at its real position p_s; the straight line
       p_s -> x* is the unfolded path, whose first segment equals the
       physical first segment. */
    V2 xp = { x, 0.0 };
    for (int j = k - 1; j >= 0; j--)
        xp = reflect_point(xp, g_mir[chain[j]].a, g_mir[chain[j]].b);

    V2 ps = { xl, g_l0.y };
    V2 dps = { 1.0, 0.0 }; /* d p_s / d xl (real light) */

    V2 delta = vsub(xp, ps);
    double r = vlen(delta);
    if (r < EPS) return 0;
    V2 s = vscale(delta, 1.0 / r);
    /* |d theta / d xl|: angle jacobian of the unfolded straight line,
       reflections are orthogonal so the magnitude carries to the floor. */
    p->dtheta = fabs(vcross(delta, dps)) / (r * r);

    /* Fold by forward tracing from the real light along the unfolded
       direction; each mirror hit reflects the direction. */
    V2 d0 = s;          /* physical direction leaving the light */
    V2 d = d0;
    V2 o = { xl, g_l0.y };
    double r0 = -1.0;
    for (int j = 0; j < k; j++) {
        double t;
        int hit = trace_mirror(o, d, -1.0, &t);
        if (hit != chain[j]) return 0;
        if (j == 0) r0 = t;
        o = vadd(o, vscale(d, t));
        d = reflect_dir(d, g_mir[hit].a, g_mir[hit].b);
    }
    /* last segment must reach the floor at x without a mirror in the way */
    if (d.y >= -1e-12) return 0;
    double tf = -o.y / d.y;
    if (tf <= EPS) return 0;
    double xf = o.x + tf * d.x;
    if (fabs(xf - x) > 1e-7) return 0;
    double mt;
    if (trace_mirror(o, d, tf, &mt) >= 0) return 0;

    p->cL = -d0.y;
    p->cosX = -d.y;
    p->r0 = k == 0 ? tf : r0;
    p->rk = tf;
    if (p->cL <= 0 || p->cosX <= 0) return 0;
    p->G0 = p->cL / p->r0;
    p->Gk = p->cosX / p->rk;
    p->ok = 1;
    return 1;
}

/* path contribution per unit light length: D_k = integral f d l */
static double path_f(const Path* p) {
    if (!p->ok) return 0.0;
    return g_Le * p->rho * p->G0 * p->Gk;
}
static double pdf_light(void) { return 1.0 / g_L; }
static double pdf_brdf(const Path* p) {
    if (!p->ok || p->dtheta <= 0) return 0.0;
    return 0.5 * p->cosX * p->dtheta;
}
/* power heuristic (exponent 2): still a partition of unity and unbiased,
   but down-weights a strategy where the other pdf dominates, which tames
   the high-variance BRDF half on narrow-support caustic classes. */
static void mis_weights(double pL, double pB, double* wL, double* wB) {
    double a = pL * pL, b = pB * pB;
    double s = a + b;
    if (s <= 0) { *wL = 1.0; *wB = 0.0; return; }
    *wL = a / s;
    *wB = b / s;
}

/* partition-of-unity gate: composition is allowed only if weights are
   non-negative and sum to one within tolerance */
static int gate_pass(const double* w, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        if (w[i] < -1e-12) return 0;
        s += w[i];
    }
    return fabs(s - 1.0) <= 1e-9;
}

/* ---- receiver-side tracing ------------------------------------------------ */

/* Cosine-weighted flatland sample (pdf = cos/2): sin(theta) = 2u-1.
   Traces the deterministic specular chain and connects to the light. */
static int brdf_trace(double x, double u, Path* p) {
    double sinth = 2.0 * u - 1.0;
    double theta = asin(sinth);
    V2 d = { sin(theta), cos(theta) };
    V2 o = { x, 0.0 };
    int chain[MAX_K], k = 0;
    double xl;
    for (int step = 0; step <= MAX_K; step++) {
        double tl;
        int light_first = trace_light(o, d, -1.0, &xl, &tl);
        double tm;
        int m = trace_mirror(o, d, light_first ? tl : -1.0, &tm);
        if (light_first && m < 0) {
            return solve_chain(x, chain, k, xl, p);
        }
        if (m < 0 || k == MAX_K) return 0;
        chain[k++] = m;
        o = vadd(o, vscale(d, tm));
        d = reflect_dir(d, g_mir[m].a, g_mir[m].b);
    }
    return 0;
}

/* ---- reference and analytic direct light --------------------------------- */

static double class_reference(double x, const int* chain, int k, int strata) {
    double acc = 0.0;
    for (int i = 0; i < strata; i++) {
        double u = ((double)i + 0.5) / strata;
        double xl = g_l0.x + g_L * u;
        Path p;
        if (solve_chain(x, chain, k, xl, &p)) acc += path_f(&p);
    }
    return g_L * acc / strata;
}

/* closed-form direct (k=0) integral, per unit Le:
   integral y^2 / ((u-x)^2+y^2)^2 du */
static double direct_analytic(double x) {
    double y = g_l0.y;
    double F = 0.0;
    double va = g_l0.x - x, vb = g_l1.x - x;
    /* antiderivative, via v=y tan t:
       F(v) = (atan(v/y) + v y/(v^2+y^2)) / (2y) */
    double Fa = (atan(va / y) + va * y / (va * va + y * y)) / (2.0 * y);
    double Fb = (atan(vb / y) + vb * y / (vb * vb + y * y)) / (2.0 * y);
    F = Fb - Fa;
    return g_Le * F;
}

/* ---- RNG ------------------------------------------------------------------ */

static unsigned long long rng_state;
static double randu(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 33) / (double)(1ULL << 31);
}

/* ---- estimators ----------------------------------------------------------- */

typedef struct { double mean, var; } Est;

static Est est_light_iid(double x, const int* chain, int k, int n) {
    Est e = { 0, 0 };
    double pL = pdf_light(), sum = 0, sum2 = 0;
    for (int i = 0; i < n; i++) {
        double xl = g_l0.x + g_L * randu();
        Path p;
        double v = 0;
        if (solve_chain(x, chain, k, xl, &p)) v = path_f(&p) / pL;
        sum += v; sum2 += v * v;
    }
    e.mean = sum / n;
    double var = n > 1 ? (sum2 / n - (sum / n) * (sum / n)) * n / (n - 1) : 0;
    e.var = var > 0 ? var : 0;
    return e;
}

static Est est_combined(double x, const int* chain, int k, int n) {
    Est e = { 0, 0 };
    double pL = pdf_light();
    double sumL = 0, sum2L = 0, sumB = 0, sum2B = 0;
    for (int i = 0; i < n; i++) {
        double xl = g_l0.x + g_L * ((i + randu()) / n);
        Path p;
        double vL = 0;
        if (solve_chain(x, chain, k, xl, &p)) {
            double pB = pdf_brdf(&p), wL, wB;
            mis_weights(pL, pB, &wL, &wB);
            if (!gate_pass((double[]){ wL, wB }, 2)) { e.mean = -1.0; return e; }
            vL = wL * path_f(&p) / pL;
        }
        sumL += vL; sum2L += vL * vL;
    }
    for (int i = 0; i < n; i++) {
        Path p;
        double vB = 0;
        if (brdf_trace(x, randu(), &p) &&
            memcmp(p.chain, chain, sizeof(int) * k) == 0 && p.k == k) {
            double pB = pdf_brdf(&p), wL, wB;
            mis_weights(pL, pB, &wL, &wB);
            if (!gate_pass((double[]){ wL, wB }, 2)) { e.mean = -1.0; return e; }
            vB = wB * path_f(&p) / pB;
        }
        sumB += vB; sum2B += vB * vB;
    }
    double mean = (sumL + sumB) / n;
    double mL = sumL / n, mB = sumB / n;
    double vL = n > 1 ? (sum2L / n - mL * mL) * n / (n - 1) : 0;
    double vB = n > 1 ? (sum2B / n - mB * mB) * n / (n - 1) : 0;
    e.mean = mean;
    e.var = (vL > 0 ? vL : 0) + (vB > 0 ? vB : 0); /* var of the sum of means */
    return e;
}

/* naive unidirectional PT: cosine ray, no connection across the specular
   chain. In flatland it occasionally lands on the light by chance; most
   caustic pixels get zero hits at practical sample counts. */
static Est est_naive(double x, const int* chain, int k, int n) {
    Est e = { 0, 0 };
    double sum = 0, sum2 = 0;
    for (int i = 0; i < n; i++) {
        Path p;
        double v = 0;
        if (brdf_trace(x, randu(), &p) &&
            memcmp(p.chain, chain, sizeof(int) * k) == 0 && p.k == k) {
            double pB = pdf_brdf(&p);
            if (pB > 0) v = path_f(&p) / pB;
        }
        sum += v; sum2 += v * v;
    }
    e.mean = sum / n;
    double var = n > 1 ? (sum2 / n - e.mean * e.mean) * n / (n - 1) : 0;
    e.var = var > 0 ? var : 0;
    return e;
}

/* ---- checks --------------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    scene_build();
    rng_state = 0x123456789abcdef0ULL;

    /* -- vector/ reflection basics -- */
    V2 d = { 1, -1 };
    V2 up = reflect_dir(d, (V2){ 0, 0 }, (V2){ 1, 0 }); /* floor line */
    V2 expect = { 1, 1 };
    check("reflect_dir", fabs(up.x - expect.x) < 1e-12 && fabs(up.y - expect.y) < 1e-12);
    V2 pr = reflect_point((V2){ 1, 2 }, (V2){ 0, 0 }, (V2){ 1, 0 });
    check("reflect_point", fabs(pr.x - 1) < 1e-12 && fabs(pr.y + 2) < 1e-12);

    double t, u;
    check("ray_seg_hit", ray_seg((V2){ 0, 1 }, (V2){ 0, -1 }, (V2){ -1, 0 }, (V2){ 1, 0 }, &t, &u)
                        && fabs(t - 1) < 1e-12 && fabs(u - 0.5) < 1e-12);
    check("ray_seg_miss", !ray_seg((V2){ 5, 1 }, (V2){ 0, -1 }, (V2){ -1, 0 }, (V2){ 1, 0 }, &t, &u));

    /* -- minimum object separation drives the deterministic Hoeffding bound -- */
    double rmin = 1e30;
    for (int i = 0; i < g_nmir; i++) {
        rmin = fmin(rmin, fmin(g_mir[i].a.y, g_mir[i].b.y)); /* mirror to floor */
        rmin = fmin(rmin, seg_seg_dist(g_l0, g_l1, g_mir[i].a, g_mir[i].b));
        for (int j = i + 1; j < g_nmir; j++)
            rmin = fmin(rmin, seg_seg_dist(g_mir[i].a, g_mir[i].b, g_mir[j].a, g_mir[j].b));
    }
    check("scene_separated", rmin > 0.15);

    /* -- enumerate mirror chains up to length MAX_K-1 -- */
    int chains[MAX_CHAIN][MAX_K];
    int clen[MAX_CHAIN], nchain = 0;
    chains[nchain][0] = 0; clen[nchain] = 0; nchain++;
    for (int k = 1; k <= MAX_K - 1; k++) {
        int total = 1;
        for (int q = 0; q < k; q++) total *= g_nmir;
        for (int code = 0; code < total; code++) {
            int c = code;
            for (int q = 0; q < k; q++) { chains[nchain][q] = c % g_nmir; c /= g_nmir; }
            clen[nchain] = k;
            if (++nchain >= MAX_CHAIN) break;
        }
    }

    /* -- unfolding round trip on a found caustic path -- */
    int roundtrip = 0, any_caustic = 0, two_bounce = 0;
    for (int ci = 0; ci < nchain; ci++) {
        if (clen[ci] == 0) continue;
        for (int pi = 0; pi < PIXELS; pi++) {
            double x = XMIN + (XMAX - XMIN) * pi / (PIXELS - 1);
            for (int i = 0; i < 64; i++) {
                double xl = g_l0.x + g_L * (i + 0.5) / 64;
                Path p;
                if (solve_chain(x, chains[ci], clen[ci], xl, &p)) {
                    roundtrip = 1;
                    any_caustic = 1;
                    if (clen[ci] >= 2) two_bounce = 1;
                    break;
                }
            }
        }
    }
    check("unfold_roundtrip", roundtrip);
    check("caustic_class_exists", any_caustic);
    check("two_bounce_exists", two_bounce);

    /* -- length-1 path degenerates to direct light -- */
    int none0[1] = { 0 };
    double ref0_q = class_reference(0.0, none0, 0, 4096);
    double ref0_a = direct_analytic(0.0);
    check("len1_direct_analytic", fabs(ref0_q - ref0_a) < 2e-4 * (1.0 + ref0_a));
    Est e0 = est_combined(0.0, none0, 0, 1024);
    check("len1_estimate_matches", fabs(e0.mean - ref0_a) < 0.05 * (1.0 + ref0_a));

    /* -- partition of unity gate -- */
    double clean[2] = { 0.37, 0.63 };
    check("gate_clean", gate_pass(clean, 2));
    double corrupt[2] = { 0.5, 0.25 };
    check("gate_blocks_corrupt", !gate_pass(corrupt, 2));
    double neg[2] = { -0.2, 1.2 };
    check("gate_blocks_negative", !gate_pass(neg, 2));

    /* -- weight partition across every accepted sample of the whole image -- */
    int wbad = 0, wcount = 0;
    for (int pi = 0; pi < PIXELS; pi++) {
        double x = XMIN + (XMAX - XMIN) * pi / (PIXELS - 1);
        for (int ci = 0; ci < nchain; ci++) {
            for (int i = 0; i < 32; i++) {
                double xl = g_l0.x + g_L * (i + 0.5) / 32;
                Path p;
                if (solve_chain(x, chains[ci], clen[ci], xl, &p)) {
                    double wL, wB;
                    mis_weights(pdf_light(), pdf_brdf(&p), &wL, &wB);
                    if (!gate_pass((double[]){ wL, wB }, 2)) wbad++;
                    wcount++;
                }
            }
        }
    }
    check("weights_partition", wcount > 100 && wbad == 0);

    /* -- reference image + error-bound coverage -- */
    const int N = 256, N4 = 1024, Q = 2048;
    const double alpha = 0.05;
    /* Deterministic upper bound on every light-strategy term f/pL:
       f = Le rho (cL/r0)(cosX/rk), pL = 1/L, cosines <= 1, every segment
       length >= rmin  =>  f/pL <= Le L / rmin^2. */
    double bound_b = g_Le * g_L / (rmin * rmin);
    /* constructive check: no observed term exceeds the declared bound */
    double term_max = 0.0; int term_n = 0;
    for (int pi = 0; pi < PIXELS; pi++) {
        double x = XMIN + (XMAX - XMIN) * pi / (PIXELS - 1);
        for (int ci = 0; ci < nchain; ci++) {
            for (int i = 0; i < 64; i++) {
                double xl = g_l0.x + g_L * randu();
                Path p;
                if (solve_chain(x, chains[ci], clen[ci], xl, &p)) {
                    double v = path_f(&p) / pdf_light();
                    if (v > term_max) term_max = v;
                    term_n++;
                }
            }
        }
    }
    check("term_under_bound", term_n > 100 && term_max <= bound_b * (1.0 + 1e-9));
    double eps_hoeff = bound_b * sqrt(log(2.0 / alpha) / (2.0 * N));
    int cells = 0, hcov = 0, ccov = 0, ccov4 = 0;
    double se_l = 0, se_b = 0, se_m = 0, se_m4 = 0, se_n = 0;
    double naive_zero = 0, caustic_cells = 0;
    double ref_img[PIXELS] = { 0 }, mis_img[PIXELS] = { 0 };

    for (int pi = 0; pi < PIXELS; pi++) {
        double x = XMIN + (XMAX - XMIN) * pi / (PIXELS - 1);
        for (int ci = 0; ci < nchain; ci++) {
            double ref = class_reference(x, chains[ci], clen[ci], Q);
            if (clen[ci] == 0) ref_img[pi] += ref;
            if (ref < 1e-6) continue;
            cells++;
            if (clen[ci] > 0) caustic_cells++;

            rng_state = 0xabc00000ULL ^ (unsigned long long)(pi * 131 + ci * 17);
            Est l = est_light_iid(x, chains[ci], clen[ci], N);
            Est m = est_combined(x, chains[ci], clen[ci], N);
            Est n4 = est_combined(x, chains[ci], clen[ci], N4);
            Est nv = est_naive(x, chains[ci], clen[ci], 8);
            if (m.mean < 0) { check("gate_inside_estimator", 0); return 1; }

            if (fabs(l.mean - ref) <= eps_hoeff) hcov++;
            double eps_cheb = sqrt(m.var / (N * alpha));
            double eps_cheb4 = sqrt(n4.var / (N4 * alpha));
            if (fabs(m.mean - ref) <= eps_cheb) ccov++;
            if (fabs(n4.mean - ref) <= eps_cheb4) ccov4++;

            se_l += (l.mean - ref) * (l.mean - ref);
            se_m += (m.mean - ref) * (m.mean - ref);
            se_m4 += (n4.mean - ref) * (n4.mean - ref);
            se_n += (nv.mean - ref) * (nv.mean - ref);
            if (clen[ci] > 0 && nv.mean == 0.0) naive_zero++;
            if (clen[ci] == 0) mis_img[pi] += m.mean;
        }
        /* add caustic class estimates into the MIS image strip */
        for (int ci = 0; ci < nchain; ci++) {
            if (clen[ci] == 0) continue;
            Est m = est_combined(x, chains[ci], clen[ci], N);
            if (m.mean > 0) mis_img[pi] += m.mean;
        }
    }

    /* BRDF-only RMSE for the variance comparison */
    for (int pi = 0; pi < PIXELS; pi++) {
        double x = XMIN + (XMAX - XMIN) * pi / (PIXELS - 1);
        for (int ci = 0; ci < nchain; ci++) {
            double ref = class_reference(x, chains[ci], clen[ci], Q);
            if (ref < 1e-6) continue;
            rng_state = 0xdef00000ULL ^ (unsigned long long)(pi * 977 + ci * 31);
            double sum = 0;
            for (int i = 0; i < N; i++) {
                Path p;
                if (brdf_trace(x, randu(), &p) &&
                    memcmp(p.chain, chains[ci], sizeof(int) * clen[ci]) == 0 &&
                    p.k == clen[ci]) {
                    double pB = pdf_brdf(&p);
                    if (pB > 0) sum += path_f(&p) / pB;
                }
            }
            double mean = sum / N;
            se_b += (mean - ref) * (mean - ref);
        }
    }

    check("hoeffding_coverage", cells > 50 && (double)hcov / cells >= 0.92);
    check("chebyshev_coverage", (double)ccov / cells >= 0.90);
    check("chebyshev_tightens_with_n", (double)ccov4 / cells >= (double)ccov / cells - 0.02);

    double rmse_l = sqrt(se_l / cells), rmse_b = sqrt(se_b / cells);
    double rmse_m = sqrt(se_m / cells), rmse_m4 = sqrt(se_m4 / cells);
    double rmse_n = sqrt(se_n / caustic_cells);
    check("mis_beats_light", rmse_m <= rmse_l * 1.05);
    check("mis_beats_brdf", rmse_m <= rmse_b * 1.05);
    check("mc_converges", rmse_m4 < rmse_m * 0.75);
    check("naive_pt_sparse", caustic_cells > 20 &&
          naive_zero / caustic_cells >= 0.5 && rmse_n > 2.0 * rmse_m);

    /* -- sanity: reference caustic energy is visible where naive PT is blind -- */
    double caus_total = 0;
    for (int ci = 0; ci < nchain; ci++) {
        if (clen[ci] == 0) continue;
        for (int pi = 0; pi < PIXELS; pi++) {
            double x = XMIN + (XMAX - XMIN) * pi / (PIXELS - 1);
            caus_total += class_reference(x, chains[ci], clen[ci], 512);
        }
    }
    check("reference_has_caustics", caus_total > 1.0);

    /* -- compact image strips: reference vs MIS composite -- */
    const char* ramp = " .:-=+*#%@";
    char strip_ref[PIXELS + 1], strip_mis[PIXELS + 1];
    double vmax = 0;
    for (int pi = 0; pi < PIXELS; pi++) {
        double r = 0;
        double x = XMIN + (XMAX - XMIN) * pi / (PIXELS - 1);
        for (int ci = 0; ci < nchain; ci++)
            r += class_reference(x, chains[ci], clen[ci], 256);
        ref_img[pi] = r;
        vmax = fmax(vmax, r);
    }
    for (int pi = 0; pi < PIXELS; pi++) {
        int ir = (int)fmin(9.0, 9.0 * ref_img[pi] / vmax);
        int im = (int)fmin(9.0, 9.0 * mis_img[pi] / vmax);
        strip_ref[pi] = ramp[ir];
        strip_mis[pi] = ramp[im];
    }
    strip_ref[PIXELS] = strip_mis[PIXELS] = 0;
    printf("ref : %s\n", strip_ref);
    printf("mis : %s\n", strip_mis);
    printf("rmin=%.4f cells=%d rmse_l=%.5f rmse_b=%.5f rmse_mis=%.5f rmse_mis4=%.5f rmse_naive=%.5f\n",
           rmin, cells, rmse_l, rmse_b, rmse_m, rmse_m4, rmse_n);
    printf("naive_zero_frac=%.3f hoeff_eps=%.5f\n", naive_zero / caustic_cells, eps_hoeff);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
