
/* ao: cosine-weighted hemisphere ambient occlusion.
   - cosine-weighted hemisphere samples around the shading normal
   - distance falloff: occluder hit at t within radius r weighs 1 - t/r
   - samples with non-positive cosine weight are rejected (resampled)
   - Monte Carlo convergence: RMSE ~ 1/sqrt(N), so doubling N cuts
     noise by 1/sqrt(2); Hoeffding and empirical Bernstein bounds give
     distribution-free error upper bounds versus analytic truth and a
     high-sample reference. */
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_SURF 8

typedef struct { float x, y, z; } V3;

static V3 add3(V3 a, V3 b) { V3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static V3 sub3(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static V3 scale3(V3 a, float s) { V3 r = { a.x * s, a.y * s, a.z * s }; return r; }
static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 cross3(V3 a, V3 b) {
    V3 r = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    return r;
}
static V3 norm3(V3 a) {
    float l = sqrtf(dot3(a, a));
    V3 r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}
static float axis3(V3 v, int a) { return a == 0 ? v.x : (a == 1 ? v.y : v.z); }

/* xorshift32, draws are [0, 1) */
static uint32_t xs32(uint32_t* s) {
    uint32_t x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    return *s = x;
}
static float randf(uint32_t* s) { return (float)(xs32(s) >> 8) * (1.0f / 16777216.0f); }

/* orthonormal tangent frame around normal n: wi = lx*t + ly*b + lz*n */
static void make_frame(V3 n, V3* t, V3* b) {
    V3 helper = fabsf(n.x) > 0.9f ? (V3){ 0, 1, 0 } : (V3){ 1, 0, 0 };
    *t = norm3(cross3(helper, n));
    *b = cross3(n, *t);
}

/* cosine-weighted hemisphere sample in local frame, pdf = cos/pi */
static V3 cosine_dir(float u1, float u2, V3 t, V3 b, V3 n) {
    float r = sqrtf(u1);
    float phi = 2.0f * (float)M_PI * u2;
    float lz = sqrtf(1.0f - u1);
    V3 local = { r * cosf(phi), r * sinf(phi), lz };
    return add3(add3(scale3(t, local.x), scale3(b, local.y)), scale3(n, local.z));
}

/* adversarial generator: uniform sphere, dips below the tangent plane */
static V3 fuzz_dir(float u1, float u2, V3 t, V3 b, V3 n) {
    float z = 1.0f - 2.0f * u1;
    float rr = sqrtf(fmaxf(0.0f, 1.0f - z * z));
    float phi = 2.0f * (float)M_PI * u2;
    V3 local = { rr * cosf(phi), rr * sinf(phi), z };
    return add3(add3(scale3(t, local.x), scale3(b, local.y)), scale3(n, local.z));
}

/* scene: infinite planes and AABB walls, nearest-hit occlusion rays */
typedef struct { int kind; int axis; float q; float lo[3], hi[3]; } Surf;
typedef struct { Surf s[MAX_SURF]; int n; } Scene;

static void add_plane(Scene* sc, int axis, float q) {
    Surf* s = &sc->s[sc->n++];
    s->kind = 0; s->axis = axis; s->q = q;
}
static void add_box(Scene* sc, const float lo[3], const float hi[3]) {
    Surf* s = &sc->s[sc->n++];
    s->kind = 1;
    for (int i = 0; i < 3; i++) { s->lo[i] = lo[i]; s->hi[i] = hi[i]; }
}

static int ray_plane(const Surf* s, V3 o, V3 d, float tmax, float* t) {
    float da = axis3(d, s->axis), oa = axis3(o, s->axis);
    if (fabsf(da) < 1e-9f) return 0;
    float tt = (s->q - oa) / da;
    if (tt <= 1e-6f || tt >= tmax) return 0;
    *t = tt;
    return 1;
}
static int ray_box(const Surf* s, V3 o, V3 d, float tmax, float* t) {
    float t0 = 1e-6f, t1 = tmax;
    for (int i = 0; i < 3; i++) {
        float da = axis3(d, i), oa = axis3(o, i);
        if (fabsf(da) < 1e-9f) {
            if (oa < s->lo[i] || oa > s->hi[i]) return 0;
        } else {
            float inv = 1.0f / da;
            float ta = (s->lo[i] - oa) * inv;
            float tb = (s->hi[i] - oa) * inv;
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    *t = t0;
    return 1;
}
static int nearest_hit(const Scene* sc, V3 o, V3 d, float tmax, float* t) {
    float best = tmax, hit = 0;
    for (int i = 0; i < sc->n; i++) {
        float tt;
        int ok = sc->s[i].kind == 0 ? ray_plane(&sc->s[i], o, d, tmax, &tt)
                                    : ray_box(&sc->s[i], o, d, tmax, &tt);
        if (ok && tt < best) { best = tt; hit = 1; }
    }
    if (hit) *t = best;
    return hit;
}

/* one AO probe. returns 0 with *occ filled when accepted, -1 when the
   cosine weight dot(wi, n) is non-positive (the sample is rejected).
   pdf = cos/pi cancels the cosine in the integral, so each accepted
   probe contributes its distance-falloff occlusion directly. */
static int ao_probe(const Scene* sc, V3 p, V3 n, V3 wi, float radius, float* occ) {
    float cosw = dot3(wi, n);
    if (cosw <= 0.0f) return -1;
    V3 o = add3(p, scale3(n, 1e-4f));
    float t;
    if (nearest_hit(sc, o, wi, radius, &t)) {
        float w = 1.0f - t / radius;
        *occ = w < 0.0f ? 0.0f : (w > 1.0f ? 1.0f : w);
    } else {
        *occ = 0.0f;
    }
    return 0;
}

typedef struct { double mean, var; int accepted, rejected; } AORun;

static AORun ao_estimate(const Scene* sc, V3 p, V3 n, float radius,
                         int samples, int fuzz, uint32_t seed) {
    V3 t, b;
    make_frame(n, &t, &b);
    double mean = 0.0, m2 = 0.0;
    int accepted = 0, rejected = 0, guard = 0;
    while (accepted < samples) {
        if (guard++ > samples * 100 + 100) break;
        float u1 = randf(&seed), u2 = randf(&seed);
        V3 wi = fuzz ? fuzz_dir(u1, u2, t, b, n)
                     : cosine_dir(u1, u2, t, b, n);
        float occ;
        if (ao_probe(sc, p, n, wi, radius, &occ) != 0) { rejected++; continue; }
        accepted++;
        double d = (double)occ - mean;
        mean += d / accepted;
        m2 += d * ((double)occ - mean);
    }
    AORun r;
    r.mean = mean;
    r.var = accepted > 1 ? m2 / (accepted - 1) : 0.0;
    r.accepted = accepted;
    r.rejected = rejected;
    return r;
}

/* analytic reference: infinite ceiling plane at distance d, AO radius r.
   A = (1/pi) int_{t<r} (1-t/r) cos dω = (1 - d/r)^2 for d < r, else 0. */
static double analytic_mean(double d, double r) {
    if (d >= r) return 0.0;
    double c0 = d / r;
    return (1.0 - c0) * (1.0 - c0);
}
/* single-probe variance for the same scene (closed form, see README) */
static double analytic_var(double d, double r) {
    double c0 = d / r;
    if (c0 >= 1.0) return 0.0;
    double e2 = 1.0 - 4.0 * c0 + 3.0 * c0 * c0 - 2.0 * c0 * c0 * log(c0);
    double m = (1.0 - c0) * (1.0 - c0);
    return e2 - m * m;
}

/* distribution-free error upper bounds for a mean of N [0,1] draws */
static double hoeff_eps(int N, double delta) {
    return sqrt(log(2.0 / delta) / (2.0 * N));
}
static double bernstein_eps(double sample_var, int N, double delta) {
    return sqrt(2.0 * sample_var * log(2.0 / delta) / N)
         + 7.0 * log(2.0 / delta) / (3.0 * (N - 1));
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define LEVELS 9
#define TRIALS 200

int main(void) {
    V3 up = { 0, 1, 0 };

    /* frame and sampler geometry */
    V3 t, b;
    make_frame(up, &t, &b);
    check("frame_orthonormal", fabsf(dot3(t, b)) < 1e-6f && fabsf(dot3(t, up)) < 1e-6f &&
                               fabsf(dot3(b, up)) < 1e-6f &&
                               fabsf(dot3(t, t) - 1.0f) < 1e-6f);
    V3 wi = cosine_dir(0.25f, 0.4f, t, b, up);
    check("cos_dir_unit", fabsf(dot3(wi, wi) - 1.0f) < 1e-5f);
    check("cos_dir_hemisphere", dot3(wi, up) >= 0.0f);

    /* negative / zero cosine weight must be rejected */
    Scene empty;
    empty.n = 0;
    float occ;
    V3 down = scale3(up, -1.0f);
    check("reject_negative", ao_probe(&empty, (V3){ 0, 1, 0 }, up, down, 2.0f, &occ) == -1);
    check("reject_zero", ao_probe(&empty, (V3){ 0, 1, 0 }, up, t, 2.0f, &occ) == -1);

    /* adversarial uniform-sphere stream: accounting and [0,1] output */
    AORun fz = ao_estimate(&empty, (V3){ 0, 1, 0 }, up, 2.0f, 1000, 1, 0x9e3779b9u);
    check("fuzz_accounting", fz.accepted == 1000 && fz.rejected > 0);
    check("fuzz_bounded", fz.mean >= 0.0 && fz.mean <= 1.0);
    /* uniform sphere dips below the tangent plane half the time */
    double rej_frac = (double)fz.rejected / (fz.accepted + fz.rejected);
    check("fuzz_rejects_half", rej_frac > 0.35 && rej_frac < 0.65);
    AORun clean = ao_estimate(&empty, (V3){ 0, 1, 0 }, up, 2.0f, 1000, 0, 0x9e3779b9u);
    check("clean_no_reject", clean.rejected == 0 && clean.accepted == 1000);

    /* distance falloff: plane 1 unit above point, radius 2 -> weight 0.5 */
    Scene ceil_scene;
    ceil_scene.n = 0;
    add_plane(&ceil_scene, 1, 2.0f);
    check("distance_weight_half",
          ao_probe(&ceil_scene, (V3){ 0, 1, 0 }, up, up, 2.0f, &occ) == 0 &&
          fabsf(occ - 0.5f) < 1e-3f);
    /* surface offset shrinks the usable ray: radius shorter than the gap misses */
    check("radius_boundary_miss",
          nearest_hit(&ceil_scene, (V3){ 0, 1.0001f, 0 }, up, 0.99f, &occ) == 0);
    check("miss_weight_zero",
          ao_probe(&empty, (V3){ 0, 1, 0 }, up, up, 2.0f, &occ) == 0 && occ == 0.0f);

    /* AABB occluder: hit through the slab, miss beside it */
    Scene box_scene;
    box_scene.n = 0;
    float lo[3] = { -1, 1.5f, -1 }, hi[3] = { 1, 2.5f, 1 };
    add_box(&box_scene, lo, hi);
    float tbox;
    check("box_hit", nearest_hit(&box_scene, (V3){ 0, 0.5f, 0 }, up, 2.0f, &tbox) == 1 &&
                     fabsf(tbox - 1.0f) < 1e-3f);
    V3 beside = norm3(sub3((V3){ 5, 0.6f, 0 }, (V3){ 0, 0.5f, 0 }));
    check("box_miss", nearest_hit(&box_scene, (V3){ 0, 0.5f, 0 }, beside, 8.0f, &tbox) == 0);

    /* analytic scene: ceiling plane distance d below the AO radius r */
    double d = 1.0, r = 2.0;
    Scene plane;
    plane.n = 0;
    add_plane(&plane, 1, (float)(1.0 + d));
    double truth = analytic_mean(d, r);
    check("analytic_truth_value", fabs(truth - 0.25) < 1e-12);

    AORun hiref = ao_estimate(&plane, (V3){ 0, 1, 0 }, up, (float)r, 16384, 0, 12345u);
    check("analytic_match_hi", fabs(hiref.mean - truth) < 0.02);
    check("estimate_bounded", hiref.mean >= 0.0 && hiref.mean <= 1.0);

    /* measured single-probe variance matches the closed form */
    double var_true = analytic_var(d, r);
    check("single_variance_matches", fabs(hiref.var - var_true) < 0.15 * var_true);

    /* other analytic points, including no-occlusion configurations */
    struct { double d, r; } cfgs[] = { { 1.0, 1.5 }, { 1.0, 1.0 }, { 1.0, 0.9 } };
    int cfg_ok = 1;
    for (int k = 0; k < 3; k++) {
        Scene sc;
        sc.n = 0;
        add_plane(&sc, 1, (float)(1.0 + cfgs[k].d));
        AORun q = ao_estimate(&sc, (V3){ 0, 1, 0 }, up, (float)cfgs[k].r, 4096, 0, 777u + k);
        if (fabs(q.mean - analytic_mean(cfgs[k].d, cfgs[k].r)) > 0.02) cfg_ok = 0;
    }
    check("analytic_points", cfg_ok);

    /* convergence study: TRIALS independent runs at each sample count */
    int N[LEVELS];
    double rmse[LEVELS], grand_mean[LEVELS];
    double hoeff_cov[LEVELS], bern_cov[LEVELS];
    printf("\n%6s %10s %10s %8s %9s %9s\n", "N", "mean", "rmse", "ratio", "hoeff_c", "bern_c");
    for (int lvl = 0; lvl < LEVELS; lvl++) {
        N[lvl] = 16 << lvl;
        double se = 0.0, bias_sum = 0.0, hcov = 0.0, bcov = 0.0;
        double eh = hoeff_eps(N[lvl], 0.05);
        for (int k = 0; k < TRIALS; k++) {
            AORun q = ao_estimate(&plane, (V3){ 0, 1, 0 }, up, (float)r,
                                  N[lvl], 0, (uint32_t)(0x1000 + lvl * 7919 + k * 104729));
            double err = q.mean - truth;
            se += err * err;
            bias_sum += q.mean;
            if (fabs(err) <= eh) hcov++;
            if (fabs(err) <= bernstein_eps(q.var, N[lvl], 0.05)) bcov++;
        }
        rmse[lvl] = sqrt(se / TRIALS);
        grand_mean[lvl] = bias_sum / TRIALS;
        hoeff_cov[lvl] = hcov / TRIALS;
        bern_cov[lvl] = bcov / TRIALS;
        printf("%6d %10.5f %10.6f %8s %9.3f %9.3f\n", N[lvl], grand_mean[lvl], rmse[lvl],
               lvl == 0 ? "-" : "", hoeff_cov[lvl], bern_cov[lvl]);
    }
    for (int lvl = 1; lvl < LEVELS; lvl++)
        printf("%6d %10s %10s %8.3f\n", N[lvl], "", "", rmse[lvl] / rmse[lvl - 1]);

    /* doubling N: noise must fall by ~1/sqrt(2) per step */
    int ratio_ok = 1;
    for (int lvl = 1; lvl < LEVELS; lvl++) {
        double q = rmse[lvl] / rmse[lvl - 1];
        if (!(q > 0.50 && q < 0.95)) ratio_ok = 0;
    }
    check("doubling_ratio_band", ratio_ok);

    /* log-log slope of RMSE vs N must sit near the theoretical -1/2 */
    double sx = 0, sy = 0;
    for (int lvl = 0; lvl < LEVELS; lvl++) { sx += log((double)N[lvl]); sy += log(rmse[lvl]); }
    sx /= LEVELS; sy /= LEVELS;
    double sxx = 0, sxy = 0;
    for (int lvl = 0; lvl < LEVELS; lvl++) {
        double x = log((double)N[lvl]) - sx, y = log(rmse[lvl]) - sy;
        sxx += x * x; sxy += x * y;
    }
    double slope = sxy / sxx;
    printf("loglog_slope=%.4f (theory -0.5)\n", slope);
    check("slope_half", slope > -0.70 && slope < -0.30);

    /* endpoint scaling: 256x samples -> ~16x less noise */
    double end_ratio = rmse[0] / rmse[LEVELS - 1];
    printf("endpoint_rmse_ratio=%.2f (theory %.2f)\n", end_ratio, sqrt((double)N[LEVELS - 1] / N[0]));
    check("endpoint_scaling", end_ratio > 10.0 && end_ratio < 24.0);

    /* empirical RMSE tracks the predicted standard error s/sqrt(N) */
    double pred_se = sqrt(var_true / N[LEVELS - 1]);
    printf("rmse_vs_predicted=%.3f\n", rmse[LEVELS - 1] / pred_se);
    check("predicted_standard_error",
          rmse[LEVELS - 1] > 0.6 * pred_se && rmse[LEVELS - 1] < 1.5 * pred_se);

    /* error upper bounds actually cover at least the nominal 1-delta rate */
    int cov_ok = 1;
    for (int lvl = 0; lvl < LEVELS; lvl++)
        if (hoeff_cov[lvl] < 0.93 || bern_cov[lvl] < 0.93) cov_ok = 0;
    check("bound_coverage", cov_ok);

    /* bias of the estimator stays inside its own standard error band */
    check("unbiased", fabs(grand_mean[0] - truth) < 0.01);

    /* corner scene: two perpendicular walls, compared against a
       high-sample reference (4 x 16384 draws) */
    Scene corner;
    corner.n = 0;
    add_plane(&corner, 0, 0.0f);
    add_plane(&corner, 2, 0.0f);
    V3 ncorner = norm3((V3){ 1, 0, 1 });
    double ref = 0.0;
    for (int k = 0; k < 4; k++)
        ref += ao_estimate(&corner, (V3){ 0.2f, 0.5f, 0.2f }, ncorner, 3.0f,
                           16384, 0, 9001u + k).mean;
    ref /= 4;

    AORun near = ao_estimate(&corner, (V3){ 0.1f, 0.5f, 0.1f }, ncorner, 3.0f, 8192, 0, 55u);
    AORun far = ao_estimate(&corner, (V3){ 4.0f, 0.5f, 4.0f }, ncorner, 3.0f, 8192, 0, 56u);
    printf("ao_near=%.4f ao_far=%.4f reference=%.4f\n", near.mean, far.mean, ref);
    check("corner_darker_than_open", near.mean > far.mean + 0.10);
    check("reference_bounded", ref > 0.0 && ref < 1.0);

    /* low-N trials versus the high-sample reference: Hoeffding upper bound */
    int within = 0;
    double eps_ref = hoeff_eps(64, 0.05);
    double low_mean = 0.0;
    for (int k = 0; k < TRIALS; k++) {
        AORun q = ao_estimate(&corner, (V3){ 0.2f, 0.5f, 0.2f }, ncorner, 3.0f,
                              64, 0, (uint32_t)(31337 + k));
        if (fabs(q.mean - ref) <= eps_ref) within++;
        low_mean += q.mean;
    }
    low_mean /= TRIALS;
    printf("reference_coverage=%d/%d eps=%.4f low_mean_err=%.4f\n",
           within, TRIALS, eps_ref, fabs(low_mean - ref));
    check("reference_bound", within >= (int)(0.93 * TRIALS));
    check("reference_converges", fabs(low_mean - ref) < 0.02);

    /* same seed must reproduce the estimate bit-for-bit */
    AORun d1 = ao_estimate(&corner, (V3){ 0.2f, 0.5f, 0.2f }, ncorner, 3.0f, 256, 0, 42u);
    AORun d2 = ao_estimate(&corner, (V3){ 0.2f, 0.5f, 0.2f }, ncorner, 3.0f, 256, 0, 42u);
    check("deterministic", d1.mean == d2.mean && d1.var == d2.var);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
