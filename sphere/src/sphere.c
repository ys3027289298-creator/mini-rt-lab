/* sphere: numerically robust ray-sphere intersection.
 *
 * Ray P(t) = O + t D (|D| = 1), sphere center C radius R.
 * With L = O - C the equation is  t^2 + 2 b t + c = 0,
 * b = L.D, c = L.L - R^2, discriminant (half form) D = b^2 - c.
 *
 * Three solvers are compared:
 *   sphere_naive   - textbook float formula, -b +/- sqrt(D), hard D < 0
 *   sphere_stable  - FMA, cancellation-free roots, tangent tolerance band
 *   sphere_ref     - same stable math in double precision (reference)
 *
 * `make test` runs the self-checks and writes CSV data to the directory given
 * as argv[1] (default ./out).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#define EPS_F 1.1920928955078125e-7f   /* 2^-23 */
/* Tangent band: D in [-tau*eps*scale, +tau*eps*scale] is treated as zero,
 * scale = b^2 + |c| being the magnitude of the terms that form D. */
#define TANGENT_TAU 8.0f
/* Forward-error coefficient validated by the randomized test below. */
#define BOUND_K 64.0f
/* Tangent-case coefficient for the merged-root displacement, likewise
 * validated empirically: displacement <= BOUND_KT * sqrt(eps * scale). */
#define BOUND_KT 4.0f

typedef struct {
    float b, c, D, scale, tol;
    int tangent;
} SphereDbg;

/* ---------- reference solver: stable math in double precision ---------- */
static int sphere_ref(const double o[3], const double d[3],
                      const double ctr[3], double R, double t[2]) {
    double L[3] = { o[0] - ctr[0], o[1] - ctr[1], o[2] - ctr[2] };
    double b = L[0] * d[0] + L[1] * d[1] + L[2] * d[2];
    double ll = L[0] * L[0] + L[1] * L[1] + L[2] * L[2];
    double c = ll - R * R;
    double D = b * b - c;
    if (D < 0.0) return 0;
    double s = sqrt(D);
    double q = -(b + (b >= 0.0 ? s : -s));   /* -(b + sgn(b)*s), no cancellation */
    double r1 = q, r2 = c / q;               /* product of roots = c */
    t[0] = fmin(r1, r2);
    t[1] = fmax(r1, r2);
    return 1;
}

/* ---------- textbook float solver (the fragile baseline) ---------- */
static int sphere_naive(const float o[3], const float d[3],
                        const float ctr[3], float R, float t[2]) {
    float L[3] = { o[0] - ctr[0], o[1] - ctr[1], o[2] - ctr[2] };
    float b = L[0] * d[0] + L[1] * d[1] + L[2] * d[2];
    float ll = L[0] * L[0] + L[1] * L[1] + L[2] * L[2];
    float c = ll - R * R;
    float D = b * b - c;
    if (D < 0.0f) return 0;                   /* hard zero: tangents get lost */
    float s = sqrtf(D);
    t[0] = -b - s;                           /* catastrophic cancellation here */
    t[1] = -b + s;
    return 1;
}

/* ---------- robust float solver ----------
 * Returns: 0 = clear miss, 1 = two roots, 2 = tangent (D clamped from the
 * tolerance band, one merged root). t[] holds near/far roots (equal in the
 * tangent case). dbg may be NULL.
 */
static int sphere_stable(const float o[3], const float d[3],
                         const float ctr[3], float R, float t[2],
                         SphereDbg* dbg) {
    float L[3] = { o[0] - ctr[0], o[1] - ctr[1], o[2] - ctr[2] };
    /* FMA keeps each dot product at one rounding per term and, crucially,
     * computes L.L - R^2 and b^2 - c with a single final rounding. */
    float b = fmaf(L[0], d[0], fmaf(L[1], d[1], L[2] * d[2]));
    float ll = fmaf(L[0], L[0], fmaf(L[1], L[1], L[2] * L[2]));
    float c = ll - R * R;
    float D = fmaf(b, b, -c);                /* b*b - c in one rounding */

    float scale = b * b + fabsf(c);
    float tol = TANGENT_TAU * EPS_F * scale;
    int tangent = 0;
    if (D < 0.0f) {
        if (D >= -tol) { D = 0.0f; tangent = 1; }   /* indistinguishable from 0 */
        else { if (dbg) *dbg = (SphereDbg){b,c,D,scale,tol,0}; return 0; }
    } else if (D <= tol) {
        tangent = 1;                                /* +0 side of the band */
    }

    if (D == 0.0f) {
        t[0] = t[1] = -b;
        if (dbg) *dbg = (SphereDbg){b,c,0.0f,scale,tol,1};
        return 2;
    }
    float s = sqrtf(D);
    /* q = -(b + sgn(b)*s): the magnitude-addition branch, always exact-sign;
     * the other root comes from c/q (Vieta), never from subtracting near-equal
     * numbers. */
    float q = -(b + copysignf(s, b));
    if (q == 0.0f) { t[0] = t[1] = -b; tangent = 1; }
    else {
        float r1 = q, r2 = c / q;
        t[0] = fminf(r1, r2);
        t[1] = fmaxf(r1, r2);
    }
    if (dbg) *dbg = (SphereDbg){b,c,D,scale,tol,tangent};
    return tangent ? 2 : 1;
}

/* ---------------- test infrastructure ---------------- */
static uint32_t rng_state = 0x12345678u;
static float randf(float lo, float hi) {
    uint32_t x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_state = x;
    return lo + (hi - lo) * ((x >> 8) & 0xffffff) / (float)0x1000000;
}

static float abserr(float got, double ref) {
    return (float)fabs((double)got - ref);
}
/* relative error normalized against max(|ref|, eps*S): stays finite when the
 * true root is near zero. */
static float relerr(float got, double ref, float S) {
    double denom = fmax(fabs(ref), (double)EPS_F * S);
    return (float)(fabs((double)got - ref) / denom);
}
/* A priori forward-error bound for the two-root case.
 * The discriminant is the difference of terms of size ~scale, so its absolute
 * error is O(eps*scale) and the root error picks up 1/sqrt(D); the rest is
 * ordinary O(eps) rounding of the roots themselves. */
static float root_bound(const SphereDbg* g, double tref, float S) {
    return BOUND_K * EPS_F *
           (g->scale / sqrtf(fmaxf(g->D, EPS_F * g->scale)) +
            (float)fabs(tref) + S);
}

typedef struct {
    long lost_hits, false_hits, two_root, tangent, clear_miss;
    long naive_lost;    /* baseline solver's lost hits, counted separately */
    float max_ratio;      /* measured err / analytic bound expression, K stripped */
    float max_clamp_ratio;
    float max_poly_res;   /* max |t^2+2bt+c| / (eps*scale) over tangent reports */
} Counters;

static void counters_reset(Counters* z) {
    Counters z0 = {0};
    *z = z0;
}

/* check one float result against the double reference */
static int check_case(int hs, const float ts[2], int hr, const double tr[2],
                      const float of[3],
                      const float ctr[3], float R, const SphereDbg* g,
                      Counters* cb) {
    int bad = 0;
    if (hr && !hs) { cb->lost_hits++; return 1; }
    if (!hr && hs) {
        /* a hit with no reference hit is legal only inside the tangent band */
        cb->false_hits++;
        if (!g->tangent) bad = 1;
        return bad;
    }
    if (!hr && !hs) { cb->clear_miss++; return 0; }

    float S = sqrtf((of[0]-ctr[0])*(of[0]-ctr[0]) +
                    (of[1]-ctr[1])*(of[1]-ctr[1]) +
                    (of[2]-ctr[2])*(of[2]-ctr[2])) + R;
    if (hs == 2) {
        cb->tangent++;
        if (hr) {
            /* merged root must sit between the two true roots and within
             * sqrt(eps*scale) of either */
            float disp = fmaxf(abserr(ts[0], tr[0]), abserr(ts[1], tr[1]));
            float unit = sqrtf(EPS_F * g->scale) + EPS_F * S;
            float ratio = disp / unit;
            if (ratio > cb->max_clamp_ratio) cb->max_clamp_ratio = ratio;
            if (disp > BOUND_KT * unit) bad = 1;
        }
        /* polynomial residual at the reported root: FMA-evaluated, this is
         * exactly the D the solver saw, hence <= the band width. It is the
         * scale-correct guarantee for a tangent report (radial distance alone
         * is not, because the band can be geometrically wide when scale is). */
        float res = fabsf(fmaf(fmaf(ts[0], ts[0], 2.0f * g->b * ts[0]),
                               1.0f, g->c));
        float rr = res / (EPS_F * g->scale);
        if (rr > cb->max_poly_res) cb->max_poly_res = rr;
        if (rr > 2.0f * TANGENT_TAU) bad = 1;
    } else {
        cb->two_root++;
        for (int i = 0; i < 2; i++) {
            float err = abserr(ts[i], tr[i]);
            float bnd = root_bound(g, tr[i], S);
            float ratio = err / (bnd / BOUND_K);
            if (ratio > cb->max_ratio) cb->max_ratio = ratio;
            if (err > bnd) bad = 1;
        }
    }
    return bad;
}

int main(int argc, char** argv) {
    const char* outdir = argc > 1 ? argv[1] : "out";
    char path[512];
    FILE* f;
    int fail = 0;

    /* ensure the output directory exists instead of crashing on a NULL FILE* */
    if (mkdir(outdir, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "cannot create output directory %s: ", outdir);
        perror("");
        return 2;
    }
    #define OPEN_CSV(name) do { \
        snprintf(path, sizeof path, "%s/" name, outdir); \
        f = fopen(path, "w"); \
        if (!f) { perror(path); return 2; } \
    } while (0)

    /* ---- Test 1: magnitude sweeps, two geometries --------------------
     * Sphere R=1 at origin, ray from (0, h, s) toward the origin.
     *   grazing (h=0.5): discriminant stays O(1) while b^2, c grow like s^2
     *     -> root error grows like eps*s^2 until the tangent band takes over.
     *   central (h=0):   discriminant ~ s^2, well conditioned
     *     -> root error grows only like eps*s.
     * The grid runs to 2^56; past the double reference's own resolution the
     * rows are marked status=3 and excluded from bound accounting. */
    Counters c1; counters_reset(&c1);
    long ref_unresolved = 0;
    OPEN_CSV("sweep.csv");
    fprintf(f, "geom,s,t_ref,abs_naive,rel_naive,abs_stable,rel_stable,bound,status\n");
    for (int gi = 0; gi < 2; gi++) {
        float h = gi == 0 ? 0.5f : 0.0f;
        const char* gname = gi == 0 ? "grazing" : "central";
        for (int k = 0; k <= 144; k++) {
            int e = -16 + k / 2;                    /* integer exponent */
            float s = ldexpf(1.0f, e) * (k % 2 ? 1.414213562f : 1.0f);  /* half-step grid */
            float ctr[3] = {0,0,0}, R = 1.0f;
            float of[3] = {0.0f, h, s};
            float inv = 1.0f / sqrtf(h*h + s*s);
            float df[3] = {0.0f, -h * inv, -s * inv};
            double od[3] = {0, h, s}, dd[3] = {0, df[1], df[2]};
            double cd[3] = {0,0,0};

            double tr[2] = {0,0}; int hr = sphere_ref(od, dd, cd, 1.0, tr);
            float tn[2] = {0,0}, ts[2] = {0,0};
            int hn = sphere_naive(of, df, ctr, R, tn);
            SphereDbg g;
            int hs = sphere_stable(of, df, ctr, R, ts, &g);
            if (!hr) {
                /* reference itself cannot resolve the geometry; mark the row
                 * but require the float solver to stay inside its band */
                ref_unresolved++;
                if (hs && !g.tangent) { fprintf(stderr, "sweep: false hit at s=%g\n", s); fail = 1; }
                fprintf(f, "%s,%.9g,nan,nan,nan,nan,nan,nan,3\n", gname, s);
                continue;
            }

            float S = sqrtf(h*h + s*s) + R;
            float bnd = root_bound(&g, tr[0], S);
            float an = hn ? abserr(tn[0], tr[0]) : NAN;
            float rn = hn ? relerr(tn[0], tr[0], S) : NAN;
            float as_ = abserr(ts[0], tr[0]);
            fprintf(f, "%s,%.9g,%.17g,%.9g,%.9g,%.9g,%.9g,%.9g,%d\n",
                    gname, s, tr[0], an, rn, as_, relerr(ts[0], tr[0], S), bnd, hs);
            if (!hn) c1.naive_lost++;
            if (check_case(hs, ts, hr, tr, of, ctr, R, &g, &c1)) {
                fprintf(stderr, "sweep: bound violated at s=%g\n", s);
                fail = 1;
            }
        }
    }
    fclose(f);
    printf("test1 sweep:      two-root=%ld tangent=%ld stable-lost=%ld naive-lost=%ld "
           "ref-unresolved=%ld maxK=%.2f\n",
           c1.two_root, c1.tangent, c1.lost_hits, c1.naive_lost, ref_unresolved, c1.max_ratio);

    /* ---- Test 2: near-tangent sweep with controllable cancellation ----
     * Ray (-B + t, h, 0), unit sphere: b = -B, c = B^2 + h^2 - 1,
     * D = 1 - h^2. With B large, D is the difference of O(B^2) terms and the
     * unresolved band widens to ~eps*B^2 in D, i.e. ~sqrt(eps)*B in root
     * distance. Both sides of D = 0 are sampled. */
    Counters c2; counters_reset(&c2);
    OPEN_CSV("tangent.csv");
    fprintf(f, "B,Drel,side,naive_miss,status,t_ref_near,t_ref_far,abs_near,abs_far,clamp_disp,poly_res\n");
    float ctr[3] = {0,0,0}, R = 1.0f;
    float Bs[3] = {1.0f, 1e3f, 1e6f};
    for (int bi = 0; bi < 3; bi++) {
        float B = Bs[bi];
        /* |D| from 1e-9 to 1, eight points per decade; D <= R^2 = 1 keeps
         * h^2 = R^2 - D non-negative without any clamping. */
        for (int k = 0; k <= 72; k++) {
            float Drel = powf(10.0f, -9.0f + k / 8.0f);
            for (int side = 0; side < 2; side++) {
                float Drealf = (side == 0 ? 1.0f : -1.0f) * Drel;
                float h2 = R * R - Drealf;
                float h = sqrtf(h2);
                float of[3] = {-B, h, 0.0f};
                float df[3] = {1.0f, 0.0f, 0.0f};
                double od[3] = {-B, (double)h, 0}, dd[3] = {1, 0, 0}, cd[3] = {0,0,0};

                double tr[2] = {0,0}; int hr = sphere_ref(od, dd, cd, 1.0, tr);
                float tn[2] = {0,0}, ts[2] = {0,0};
                int hn = sphere_naive(of, df, ctr, R, tn);
                SphereDbg g;
                int hs = sphere_stable(of, df, ctr, R, ts, &g);

                float an = NAN, af = NAN, disp = 0.0f, polyres = NAN;
                if (hr) {
                    an = abserr(ts[0], tr[0]);
                    af = abserr(ts[1], tr[1]);
                    if (hs == 2) disp = fmaxf(an, af);
                }
                if (hs) {
                    polyres = fabsf(fmaf(fmaf(ts[0], ts[0], 2.0f * g.b * ts[0]),
                                         1.0f, g.c));
                }
                fprintf(f, "%.9g,%.9g,%d,%d,%d,%.17g,%.17g,%.9g,%.9g,%.9g,%.9g\n",
                        B, Drel, side, hn ? 0 : 1, hs,
                        hr ? tr[0] : NAN, hr ? tr[1] : NAN, an, af, disp, polyres);
                if (hr && !hn) c2.naive_lost++;
                if (check_case(hs, ts, hr, tr, of, ctr, R, &g, &c2)) {
                    fprintf(stderr, "tangent: bound violated B=%g Drel=%g side=%d\n",
                            B, Drel, side);
                    fail = 1;
                }
            }
        }
    }
    fclose(f);
    printf("test2 tangent:    two-root=%ld tangent=%ld clear-miss=%ld stable-lost=%ld naive-lost=%ld false=%ld\n"
           "                  max clamp ratio=%.2f max poly-residual/(eps*scale)=%.2f\n",
           c2.two_root, c2.tangent, c2.clear_miss, c2.lost_hits, c2.naive_lost, c2.false_hits,
           c2.max_clamp_ratio, c2.max_poly_res);

    /* ---- Test 3: randomized rays over a wide magnitude range ---------- */
    Counters c3; counters_reset(&c3);
    long ref_hits = 0, ref_misses = 0, naive_lost = 0;
    float max_naive_ratio = 0.0f;
    OPEN_CSV("random.csv");
    fprintf(f, "mag,t_ref,rel_naive,rel_stable,naive_miss,status\n");
    for (long i = 0; i < 200000; i++) {
        float mag = ldexpf(1.0f, (int)randf(-4.0f, 20.0f));
        float c[3] = { randf(-1,1)*mag, randf(-1,1)*mag, randf(-1,1)*mag };
        float R = mag * randf(0.1f, 10.0f);
        float of[3] = { randf(-1,1)*mag, randf(-1,1)*mag, randf(-1,1)*mag };
        float z = randf(-1.0f, 1.0f), phi = randf(0.0f, 6.2831853f);
        float rr = sqrtf(fmaxf(0.0f, 1.0f - z*z));
        float df[3] = { rr * cosf(phi), rr * sinf(phi), z };

        double od[3] = {of[0],of[1],of[2]}, dd[3] = {df[0],df[1],df[2]};
        double cd[3] = {c[0],c[1],c[2]};
        double tr[2] = {0,0}; int hr = sphere_ref(od, dd, cd, (double)R, tr);
        float tn[2] = {0,0}, ts[2] = {0,0};
        int hn = sphere_naive(of, df, c, R, tn);
        SphereDbg g;
        int hs = sphere_stable(of, df, c, R, ts, &g);

        if (!hr) { ref_misses++; continue; }
        ref_hits++;
        float Lv[3] = {of[0]-c[0], of[1]-c[1], of[2]-c[2]};
        float S = sqrtf(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]) + R;
        float rn = hn ? relerr(tn[0], tr[0], S) : NAN;
        float rs = hs ? relerr(ts[0], tr[0], S) : NAN;
        fprintf(f, "%.9g,%.17g,%.9g,%.9g,%d,%d\n",
                mag, tr[0], rn, rs, hn ? 0 : 1, hs);
        if (!hn) naive_lost++;
        else {
            float ratio = abserr(tn[0], tr[0]) /
                          (EPS_F * (g.scale / sqrtf(fmaxf(g.D, EPS_F*g.scale)) +
                                    (float)fabs(tr[0]) + S));
            if (ratio > max_naive_ratio) max_naive_ratio = ratio;
        }
        if (check_case(hs, ts, hr, tr, of, c, R, &g, &c3)) {
            fprintf(stderr, "random: bound violated i=%ld\n", i);
            fail = 1;
        }
    }
    fclose(f);
    printf("test3 random:     ref hits=%ld misses=%ld naive-lost=%ld (naive maxK=%.1f)\n"
           "                  stable two-root=%ld tangent=%ld lost=%ld false=%d maxK=%.2f\n",
           ref_hits, ref_misses, naive_lost, max_naive_ratio,
           c3.two_root, c3.tangent, c3.lost_hits, (int)c3.false_hits, c3.max_ratio);

    /* ---- Test 4: rays originating inside the sphere (c < 0) ----------- */
    long inside_fail = 0;
    for (long i = 0; i < 100000; i++) {
        float R = randf(0.5f, 5.0f);
        float rad = R * sqrtf(randf(0.0f, 0.9f));
        float z = randf(-1,1), phi = randf(0, 6.2831853f);
        float rr = sqrtf(fmaxf(0.0f, 1-z*z));
        float of[3] = { rad*rr*cosf(phi), rad*rr*sinf(phi), rad*z };
        float z2 = randf(-1,1), phi2 = randf(0, 6.2831853f);
        float rr2 = sqrtf(fmaxf(0.0f, 1-z2*z2));
        float df[3] = { rr2*cosf(phi2), rr2*sinf(phi2), z2 };
        float ctr4[3] = {0,0,0};
        double od[3]={of[0],of[1],of[2]}, dd[3]={df[0],df[1],df[2]}, cd[3]={0,0,0};
        double tr[2] = {0,0}; int hr = sphere_ref(od, dd, cd, (double)R, tr);
        SphereDbg g; float ts[2] = {0,0};
        int hs = sphere_stable(of, df, ctr4, R, ts, &g);
        if (!hr || !hs || tr[1] <= 0.0 || ts[1] <= 0.0f ||
            relerr(ts[1], tr[1], R) > 1e-5f)
            inside_fail++;
    }
    printf("test4 inside:     failures=%ld\n", inside_fail);

    /* ---- Test 5: analytic exact tangents across magnitudes ------------
     * Ray (-m+t, m, 0) vs sphere R=m: exact tangent at t=m, D=0. */
    long tangent_exact_fail = 0;
    for (int k = -2; k <= 18; k++) {
        float m = ldexpf(1.0f, k);
        float ctr5[3] = {0,0,0}, R = m;
        float of[3] = { -m, m, 0.0f }, df[3] = { 1.0f, 0.0f, 0.0f };
        SphereDbg g; float ts[2] = {0,0};
        int hs = sphere_stable(of, df, ctr5, R, ts, &g);
        double od[3]={-m,m,0}, dd[3]={1,0,0}, cd[3]={0,0,0};
        double tr[2] = {0,0}; int hr = sphere_ref(od, dd, cd, (double)m, tr);
        if (!hr || hs == 0 || relerr(ts[0], tr[0], 2*m) > 1e-5f)
            tangent_exact_fail++;
    }
    printf("test5 exact-tan:  failures=%ld\n", tangent_exact_fail);

    long total_lost = c1.lost_hits + c2.lost_hits + c3.lost_hits;
    long total_naive = c1.naive_lost + c2.naive_lost + naive_lost;
    long total_false = c2.false_hits + c3.false_hits;
    if (total_lost != 0) { fprintf(stderr, "FAIL: stable solver lost %ld reference hits\n", total_lost); fail = 1; }
    if (inside_fail || tangent_exact_fail) { fprintf(stderr, "FAIL: special-case test failed\n"); fail = 1; }
    if (c1.max_ratio > BOUND_K || c2.max_clamp_ratio > BOUND_KT ||
        c3.max_ratio > BOUND_K || c2.max_poly_res > 2.0f * TANGENT_TAU) {
        fprintf(stderr, "FAIL: measured error exceeds analytic bound\n"); fail = 1;
    }
    printf("STABLE lost=%ld (must be 0); NAIVE lost=%ld; band-only false hits=%ld (all tangent)\n",
           total_lost, total_naive, total_false);
    printf("RESULT: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}
