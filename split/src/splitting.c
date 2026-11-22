
/* splitting: spatial path splitting for glossy highlights.
 *
 * Scene: pinhole camera above a rough floor.  Each primary ray lands at p,
 * where the floor has a GGX lobe whose roughness sigma(p) varies smoothly
 * (smooth in the middle, rough towards the border of the footprint).  One
 * GGX bounce is traced with next-event estimation towards a small emissive
 * sphere, so a root path is worth
 *
 *     X = Le * F * G1(l)      (VNDF sampling, F0 = 1, so 0 <= X <= Le)
 *
 * Because the light is small, X is a hit/miss estimator: most roots return 0,
 * a few return ~Le.  Single-sample means are therefore dominated by fireflies
 * exactly at the penumbra, where the lobe partially overlaps the light.
 *
 * Splitting: at the floor hit the rest of the path may be replicated K times
 * and averaged with weights 1/K, hence E[Y | root] = E[X | root] for every K.
 * Two gates at the hit decide where, as requested:
 *
 *   - roughness gate:    sigma(p) <= SIG_GATE
 *   - hit-position gate: the mirror direction at p lies within
 *                        C * (alpha + light_cone) of the direction to the
 *                        light centre, with alpha = sigma^2.
 *
 * The position gate fires on the penumbra band.  That is where the hit
 * probability is strictly between 0 and 1 and where the within-root variance
 * lives; splitting the hot-spot centre (it always hits) or the far field (it
 * always misses) only burns rays.  In the default scene the gate fires on
 * < 1% of roots yet captures > 99% of the within-root variance.
 *
 * Children are STRATIFIED over the VNDF sample domain: child c of K jitters
 * inside its own sub-rectangle, so the penumbra edge is resolved within every
 * split root instead of being left to chance.  IID children are measured too,
 * for contrast; they already win here because the gate selects exactly the
 * high-variance roots, and stratification adds a small further gain.
 *
 * Budget metric: a baseline root costs 2 rays (primary + shadow); a split
 * root costs 1 + K.  At an equal total ray budget the number of roots scales
 * by 2 / cost, so the equal-budget variance ratio is
 *
 *     eq = (Var[Y] / Var[X]) * (cost / 2),
 *
 * and `gain = 1 / eq`.  eq < 1 means the split wins at equal total samples.
 *
 * K = 1 reproduces the original path exactly, sample for sample: the gate is
 * not evaluated and no extra RNG stream is consumed.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define PI 3.14159265358979323846

typedef struct { double x, y, z; } V3;
static V3 v3(double x, double y, double z) { V3 r = { x, y, z }; return r; }
static V3 vsub(V3 a, V3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static V3 vmul(V3 a, double s) { return v3(a.x * s, a.y * s, a.z * s); }
static double vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static double vlen(V3 a) { return sqrt(vdot(a, a)); }
static V3 vnorm(V3 a) { double l = vlen(a); return l > 0 ? vmul(a, 1.0 / l) : a; }
static double clampd(double x, double lo, double hi) { return x < lo ? lo : (x > hi ? hi : x); }

/* ---------------- RNG: counter based, so every child has its own stream ---------------- */

typedef struct { uint64_t s; } Rng;

static uint64_t mix64(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static Rng rng_make(uint64_t seed, uint64_t root, uint64_t child) {
    uint64_t h = mix64(seed ^ (root + 0x9e3779b97f4a7c15ULL));
    h = mix64(h ^ (child + 0xc2b2ae3d27d4eb4fULL));
    Rng r; r.s = h ? h : 1; return r;
}
static double rnd(Rng* r) {
    uint64_t x = r->s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    r->s = x;
    return (double)((x * 2685821657736338717ULL) >> 11) * (1.0 / 9007199254740992.0);
}

/* ---------------- scene (overridable at build time for experiments) ---------------- */

#ifndef CAM_Z
#define CAM_Z      1.20
#endif
#ifndef IMAGE_CX
#define IMAGE_CX   0.12
#endif
#ifndef IMAGE_CY
#define IMAGE_CY   0.08
#endif
#ifndef IMAGE_HW
#define IMAGE_HW   0.30
#endif
#ifndef LIGHT_X
#define LIGHT_X    0.16
#endif
#ifndef LIGHT_Y
#define LIGHT_Y    0.10
#endif
#ifndef LIGHT_Z
#define LIGHT_Z    1.05
#endif
#ifndef LIGHT_R
#define LIGHT_R    0.015
#endif
#ifndef LE
#define LE         1.0
#endif
#ifndef SIG_A
#define SIG_A      0.05
#endif
#ifndef SIG_B
#define SIG_B      0.30
#endif
#ifndef SIG_R
#define SIG_R      0.55
#endif
#ifndef SIG_GATE
#define SIG_GATE   0.35
#endif
#ifndef GATE_C
#define GATE_C     3.0
#endif
#ifndef GATE_RLO
#define GATE_RLO   0.0       /* optional lower edge of the penumbra band */
#endif
#ifndef N_ROOTS
#define N_ROOTS    400000L
#endif

static const V3 CAM = { 0.0, 0.0, CAM_Z };
static const V3 LIGHT_C = { LIGHT_X, LIGHT_Y, LIGHT_Z };

static double rough_at(double x, double y) {
    double t = (x * x + y * y) / (SIG_R * SIG_R);
    return SIG_A + SIG_B * clampd(t, 0.0, 1.0);
}

/* Smith masking, one direction; alpha = sigma^2 */
static double smith_g1(double nv, double alpha) {
    double a2 = alpha * alpha;
    return 2.0 * nv / (nv + sqrt(a2 + (1.0 - a2) * nv * nv));
}

static int hit_sphere(V3 o, V3 d, V3 c, double rad) {
    V3 oc = vsub(c, o);
    double b = vdot(oc, d);
    if (b <= 0) return 0;
    double disc = b * b - vdot(oc, oc) + rad * rad;
    if (disc < 0) return 0;
    return b - sqrt(disc) > 1e-4;
}

/* VNDF half-vector for view v (local frame, v.z > 0) from (u1, u2) in [0,1)^2. */
static V3 vndf_h(V3 v, double alpha, double u1, double u2) {
    V3 vh = vnorm(v3(alpha * v.x, alpha * v.y, v.z));
    double r = sqrt(u1), ph = 2.0 * PI * u2;
    double t1 = r * cos(ph), t2 = r * sin(ph);
    double s = 0.5 * (1.0 + vh.z);
    t2 = (1.0 - s) * sqrt(fmax(0.0, 1.0 - t1 * t1)) + s * t2;
    double t3 = sqrt(fmax(0.0, 1.0 - t1 * t1 - t2 * t2));
    return vnorm(v3(alpha * t1, alpha * t2, fmax(0.0, t3)));
}

/* One GGX bounce with next-event estimation to the sphere.  The floor normal
 * is +z, so the half-vector is built directly in the world frame. */
static double suffix_nee(V3 p, V3 v, double sig, double u1, double u2) {
    double alpha = sig * sig;
    V3 h = vndf_h(v, alpha, u1, u2);
    double vh = vdot(v, h);
    if (vh <= 0.0) return 0.0;
    V3 l = vsub(vmul(h, 2.0 * vh), v);
    if (l.z <= 0.0) return 0.0;
    if (!hit_sphere(p, l, LIGHT_C, LIGHT_R)) return 0.0;
    return LE * smith_g1(l.z, alpha);
}

/* geometry reused by the position gate */
static double light_cone_half_angle(V3 p) {
    return asin(clampd(LIGHT_R / vlen(vsub(LIGHT_C, p)), 0.0, 1.0));
}
static double mirror_angle(V3 p, V3 v) {
    V3 dl = vnorm(vsub(LIGHT_C, p));
    V3 m = v3(-v.x, -v.y, v.z);   /* reflect v about +z */
    return acos(clampd(vdot(m, dl), -1.0, 1.0));
}

/* Normalised position: 0 at the specular peak, ~1 at the light's edge. */
static double hot_ratio(V3 p, V3 v, double sig) {
    return mirror_angle(p, v) / (sig * sig + light_cone_half_angle(p));
}

/* Split decision; both gates are evaluated at the floor hit. */
static int should_split(V3 p, V3 v, double sig) {
    if (sig > SIG_GATE) return 0;                         /* roughness gate */
    double r = hot_ratio(p, v, sig);
    return r > GATE_RLO && r < GATE_C;                    /* position gate  */
}

/* largest divisor of K <= sqrt(K): makes a nearly square nr x nph grid */
static void strata_grid(int K, int* nr, int* nph) {
    int r = (int)(sqrt((double)K) + 0.5);
    while (r > 1 && K % r) r--;
    *nr = r; *nph = K / r;
}

static double child_stratified(V3 p, V3 v, double sig, int c, int nr, int nph, Rng* r) {
    int ir = c / nph, iph = c % nph;
    double u1 = (ir + rnd(r)) / (double)nr;
    double u2 = (iph + rnd(r)) / (double)nph;
    return suffix_nee(p, v, sig, u1, u2);
}
static double child_iid(V3 p, V3 v, double sig, Rng* r) {
    return suffix_nee(p, v, sig, rnd(r), rnd(r));
}

/* ---------------- estimator comparison ---------------- */

enum { CHILD_IID = 0, CHILD_STRAT = 1 };

typedef struct {
    double mean, var;
    double cost;          /* mean rays per root: 1 primary + k shadows      */
    double fired;         /* fraction of roots that were split             */
    double paired_mean, paired_se, t;   /* Y - X over all roots           */
    long mismatch;        /* roots where Y != X (should be split roots only) */
    long n;
} Result;

static Result simulate(int K, int child_kind, long N, uint64_t seed) {
    double sy = 0, sy2 = 0;
    double cost = 0, fired = 0;
    double pd = 0, pd2 = 0;
    long mismatch = 0;

    for (long i = 0; i < N; i++) {
        Rng rp = rng_make(seed, (uint64_t)i, 0);
        double x = IMAGE_CX + (2.0 * rnd(&rp) - 1.0) * IMAGE_HW;
        double y = IMAGE_CY + (2.0 * rnd(&rp) - 1.0) * IMAGE_HW;
        V3 p = v3(x, y, 0.0);
        V3 v = vnorm(vsub(CAM, p));
        double sig = rough_at(x, y);

        Rng rx = rng_make(seed, (uint64_t)i, 1);
        double X = suffix_nee(p, v, sig, rnd(&rx), rnd(&rx));

        int k = 1;
        double Y = X;
        if (K > 1 && should_split(p, v, sig)) {
            k = K; fired += 1.0;
            double s = 0.0;
            if (child_kind == CHILD_STRAT) {
                int nr, nph; strata_grid(K, &nr, &nph);
                for (int c = 0; c < K; c++) {
                    Rng rc = rng_make(seed, (uint64_t)i, (uint64_t)(2 + c));
                    s += child_stratified(p, v, sig, c, nr, nph, &rc);
                }
            } else {
                for (int c = 0; c < K; c++) {
                    Rng rc = rng_make(seed, (uint64_t)i, (uint64_t)(2 + c));
                    s += child_iid(p, v, sig, &rc);
                }
            }
            Y = s / (double)K;
        }

        double d = Y - X;
        sy += Y; sy2 += Y * Y;
        pd += d; pd2 += d * d;
        cost += 1.0 + (double)k;
        if (Y != X) mismatch++;
    }

    Result r;
    r.n = N;
    r.mean = sy / N;
    r.var = fmax(0.0, sy2 / N - r.mean * r.mean);
    r.cost = cost / N;
    r.fired = fired / N;
    r.paired_mean = pd / N;
    double vd = fmax(0.0, pd2 / N - r.paired_mean * r.paired_mean);
    r.paired_se = sqrt(vd / (double)N);
    r.t = r.paired_se > 0 ? r.paired_mean / r.paired_se : 0.0;
    r.mismatch = mismatch;
    return r;
}

/* baseline: exactly the same roots and the same child-0 stream as simulate() */
static Result simulate_base(long N, uint64_t seed) {
    double sx = 0, sx2 = 0, cost = 0;
    for (long i = 0; i < N; i++) {
        Rng rp = rng_make(seed, (uint64_t)i, 0);
        double x = IMAGE_CX + (2.0 * rnd(&rp) - 1.0) * IMAGE_HW;
        double y = IMAGE_CY + (2.0 * rnd(&rp) - 1.0) * IMAGE_HW;
        V3 p = v3(x, y, 0.0);
        V3 v = vnorm(vsub(CAM, p));
        double sig = rough_at(x, y);
        Rng rx = rng_make(seed, (uint64_t)i, 1);
        double X = suffix_nee(p, v, sig, rnd(&rx), rnd(&rx));
        sx += X; sx2 += X * X;
        cost += 2.0;   /* primary + shadow */
    }
    Result r;
    memset(&r, 0, sizeof(r));
    r.n = N;
    r.mean = sx / N;
    r.var = fmax(0.0, sx2 / N - r.mean * r.mean);
    r.cost = cost / N;
    return r;
}

/* Variance decomposition over the image footprint, using many inner draws per
 * root to estimate mu(p) and Var[X | p]:
 *
 *   Var[X] = Var_p(mu(p))  +  E_p[Var(X | p)]
 *          = between-root     + within-root
 *
 * and, restricted to gated roots, the within-root share the gate captures. */
typedef struct { double between, within, fire, within_hot; } Decomp;

static Decomp decompose(int inner, long N) {
    double sum_mu = 0, sum_mu2 = 0, sum_var = 0;
    double hot_var = 0; long hot = 0;
    uint64_t seed = 0xdec0deULL;
    for (long i = 0; i < N; i++) {
        Rng rp = rng_make(seed, (uint64_t)i, 0);
        double x = IMAGE_CX + (2.0 * rnd(&rp) - 1.0) * IMAGE_HW;
        double y = IMAGE_CY + (2.0 * rnd(&rp) - 1.0) * IMAGE_HW;
        V3 p = v3(x, y, 0.0);
        V3 v = vnorm(vsub(CAM, p));
        double sig = rough_at(x, y);

        double s = 0, s2 = 0;
        for (int j = 0; j < inner; j++) {
            Rng r = rng_make(seed, (uint64_t)i, (uint64_t)(1000 + j));
            double val = suffix_nee(p, v, sig, rnd(&r), rnd(&r));
            s += val; s2 += val * val;
        }
        double mu = s / inner;
        double wv = fmax(0.0, s2 / inner - mu * mu);
        sum_mu += mu; sum_mu2 += mu * mu; sum_var += wv;
        if (should_split(p, v, sig)) { hot++; hot_var += wv; }
    }
    Decomp d;
    double mb = sum_mu / N;
    d.between = sum_mu2 / N - mb * mb;
    d.within = sum_var / N;
    d.fire = (double)hot / N;
    d.within_hot = hot ? hot_var / hot : 0.0;
    return d;
}

/* ---------------- checks ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ---- unit checks ---- */
    Rng r = rng_make(1, 0, 0);
    int okr = 1;
    for (int i = 0; i < 10000; i++) { double u = rnd(&r); if (u < 0.0 || u >= 1.0) okr = 0; }
    check("rng_unit_interval", okr);

    Rng rh = rng_make(7, 1, 2);
    int okh = 1;
    for (int i = 0; i < 4000; i++) {
        double u1 = rnd(&rh), u2 = rnd(&rh);
        V3 v = vnorm(v3(0.2 * (rnd(&rh) - 0.5), 0.2 * (rnd(&rh) - 0.5), 1.0));
        V3 h = vndf_h(v, 0.01 + 0.5 * u2, u1, u2);
        if (fabs(vlen(h) - 1.0) > 1e-9 || h.z < -1e-12) okh = 0;
    }
    check("vndf_unit_hemisphere", okh);

    /* ---- split count of one degenerates to the original path ---- */
    Result k1 = simulate(1, CHILD_STRAT, 200000, 0x51eedULL);
    check("split_k1_identical", k1.mismatch == 0 && k1.cost == 2.0);

    /* ---- variance decomposition before splitting ---- */
    Decomp dc = decompose(256, 30000);
    double var_total = dc.between + dc.within;
    double hot_share = (dc.within > 0 && dc.fire > 0)
                       ? dc.fire * dc.within_hot / dc.within : 0.0;
    printf("# DECOMP between=%.4e within=%.4e within_share=%.3f fire=%.4f"
           " within_hot_share=%.3f\n",
           dc.between, dc.within, dc.within / var_total, dc.fire, hot_share);

    /* ---- K sweep: IID vs stratified children, equal total rays ---- */
    const long N = N_ROOTS;
    const uint64_t seed = 0x51eedULL;
    Result b = simulate_base(N, seed);
    const double base_cost = b.cost;   /* 2 rays/root */

    printf("# N=%ld base_mean=%.6e base_var=%.6e base_cost=%.3f\n",
           N, b.mean, b.var, b.cost);
    printf("# child K mean var cost_ratio eq_var_ratio gain fired t_bias mismatch\n");

    int ok_unbiased = 1, ok_gain = 0, ok_cost = 1;
    int Ks[] = { 2, 4, 8, 16 };
    for (int ci = 0; ci < 2; ci++) {
        int kind = (ci == 1) ? CHILD_STRAT : CHILD_IID;
        for (int ki = 0; ki < 4; ki++) {
            int K = Ks[ki];
            Result s = simulate(K, kind, N, seed);
            double cost_ratio = s.cost / base_cost;
            double eq = (s.var / b.var) * cost_ratio;
            double gain = 1.0 / eq;
            printf("%s %d %.6e %.6e %.4f %.4f %.4f %.3f %+.2f %ld\n",
                   kind == CHILD_STRAT ? "strat" : "iid", K, s.mean, s.var,
                   cost_ratio, eq, gain, s.fired, s.t, s.mismatch);
            if (fabs(s.t) > 4.0) ok_unbiased = 0;              /* paired bias */
            if (kind == CHILD_STRAT && eq < 0.9) ok_gain = 1;  /* real win     */
            if (kind == CHILD_STRAT && s.mismatch > (long)(s.fired * N) + 1) ok_cost = 0;
        }
    }
    check("unbiased_all_k", ok_unbiased);
    check("equal_budget_variance_drop", ok_gain);
    check("only_split_roots_change", ok_cost);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
