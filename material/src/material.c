
/* material: Lambert + Cook-Torrance GGX BRDF.
   Sampling: cosine-weighted hemisphere and GGX NDF importance sampling,
   with analytic-PDF histogram / CDF checks, the BRDF*cos/PDF estimator
   identity and an equal-sample variance comparison. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct { float x, y, z; } V3;

static const float PI = 3.14159265358979f;

static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 norm3(V3 a) {
    float l = sqrtf(dot3(a, a));
    V3 r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}
static V3 add3(V3 a, V3 b) { V3 r = { a.x+b.x, a.y+b.y, a.z+b.z }; return r; }
static V3 sub3(V3 a, V3 b) { V3 r = { a.x-b.x, a.y-b.y, a.z-b.z }; return r; }
static V3 scale3(V3 a, float s) { V3 r = { a.x*s, a.y*s, a.z*s }; return r; }
/* reflect unit vector v about unit normal n: r = 2(v.n)n - v */
static V3 reflect3(V3 v, V3 n) { return sub3(scale3(n, 2.0f * dot3(v, n)), v); }
static int finite_f(float x) { return x == x && x < INFINITY && x > -INFINITY; }

/* Lambert: albedo / pi */
static float brdf_lambert(float albedo) { return albedo / PI; }

/* GGX / Trowbridge-Reitz normal distribution.
   D(h) = a^2 / (pi * ((n.h)^2 (a^2 - 1) + 1)^2),  a = rough^2.
   rough = 0 is a delta distribution D = delta(h - n), which has no
   ordinary density: return 0 away from the pole; the sampling path
   handles the delta case explicitly. */
static float ggx_ndf(float nh, float rough) {
    if (nh <= 0.0f) return 0.0f;
    float a = rough * rough;
    if (a <= 0.0f) return 0.0f;
    float a2 = a * a;
    float d = nh * nh * (a2 - 1.0f) + 1.0f;
    return a2 / (PI * d * d);
}

/* Smith-GGX masking/shadowing, Walter 2007:
   G1(v) = n.v / (n.v + sqrt(a^2 + (1-a^2)(n.v)^2)), G = G1(i) G1(o) */
static float smith_g1(float nv, float a) {
    if (nv <= 0.0f) return 0.0f;
    float a2 = a * a;
    return nv / (nv + sqrtf(a2 + (1.0f - a2) * nv * nv));
}
static float smith_g(float ni, float no, float a) {
    return smith_g1(ni, a) * smith_g1(no, a);
}

/* Schlick Fresnel, cosine term is v.h (equals l.h by construction) */
static float fresnel_schlick(float vh, float f0) {
    float t = 1.0f - fminf(fmaxf(vh, 0.0f), 1.0f);
    return f0 + (1.0f - f0) * t * t * t * t * t;
}

/* Cook-Torrance specular BRDF:
   f(wi,wo) = D(h) G(wi,wo) F(h,wo) / (4 cos_i cos_o), h = norm(wi+wo).
   Undefined as an ordinary function at rough = 0 (mirror limit). */
static float brdf_ggx(V3 wi, V3 wo, float rough, float f0) {
    float ci = wi.z, co = wo.z;
    if (ci <= 0.0f || co <= 0.0f || rough <= 0.0f) return 0.0f;
    V3 h = norm3(add3(wi, wo));
    float nh = h.z, vh = dot3(wo, h);
    if (nh <= 0.0f || vh <= 0.0f) return 0.0f;
    float a = rough * rough;
    float d = ggx_ndf(nh, rough);
    float g = smith_g(ci, co, a);
    float f = fresnel_schlick(vh, f0);
    return d * g * f / (4.0f * ci * co);
}

/* cosine-weighted hemisphere: pdf = cos/pi */
static float cosine_pdf(float cos_theta) { return cos_theta / PI; }

static V3 cosine_sample(float u1, float u2) {
    float r = sqrtf(u1);
    float phi = 2.0f * PI * u2;
    V3 v = { r * cosf(phi), r * sinf(phi), sqrtf(1.0f - u1 > 0 ? 1.0f - u1 : 0) };
    return norm3(v);
}

/* ---- GGX importance sampling ----------------------------------------------
   D is normalized w.r.t. projected solid angle:
       integral_hemisphere D(h) (n.h) d omega_h = 1
   so the half-vector pdf is p_h(h) = D(h) (n.h).
   Isotropic lobe => phi_h = 2 pi u2, and the theta marginal inverts as
       u1 = tan^2(theta) / (a^2 + tan^2(theta))
          => cos(theta) = sqrt((1-u1) / (1 + (a^2-1) u1))
   rough = 0: theta = 0 for every u1, i.e. h = n deterministically. */
static V3 ggx_sample_half(float u1, float u2, float rough) {
    V3 h = { 0.0f, 0.0f, 1.0f };
    if (rough <= 0.0f) return h;              /* delta lobe: h = n */
    float a2 = rough * rough; a2 = a2 * a2;  /* a^2 where a = rough^2 */
    float phi = 2.0f * PI * u2;
    float cos2 = (1.0f - u1) / (1.0f + (a2 - 1.0f) * u1);
    float cos_t = sqrtf(fmaxf(0.0f, cos2));
    float sin_t = sqrtf(fmaxf(0.0f, 1.0f - cos2));
    h.x = sin_t * cosf(phi);
    h.y = sin_t * sinf(phi);
    h.z = cos_t;
    return h;
}

/* pdf w.r.t. half-vector solid angle */
static float ggx_half_pdf(V3 h, float rough) {
    float nh = h.z;
    if (nh <= 0.0f || rough <= 0.0f) return 0.0f;  /* delta: no finite pdf */
    return ggx_ndf(nh, rough) * nh;
}

/* pdf w.r.t. incoming-direction solid angle.  The reflection map
   wi = 2(wo.h)h - wo has Jacobian |d omega_h / d omega_i| = 1/(4 wo.h),
   so p(wi) = D(h)(n.h) / (4 wo.h). */
static float ggx_dir_pdf(V3 wi, V3 wo, float rough) {
    if (wi.z <= 0.0f || wo.z <= 0.0f || rough <= 0.0f) return 0.0f;
    V3 h = norm3(add3(wi, wo));
    float nh = h.z, vh = dot3(wo, h);
    if (nh <= 0.0f || vh <= 0.0f) return 0.0f;
    return ggx_ndf(nh, rough) * nh / (4.0f * vh);
}

/* sample wi for view direction wo; *pdf = 0 marks a rejected sample
   (back-facing wi).  rough = 0 takes the explicit mirror branch. */
static V3 ggx_sample_wi(V3 wo, float u1, float u2, float rough, float *pdf) {
    V3 n = { 0.0f, 0.0f, 1.0f };
    *pdf = 0.0f;
    if (rough <= 0.0f) return reflect3(wo, n);   /* delta: pdf not finite */
    V3 h = ggx_sample_half(u1, u2, rough);
    float vh = dot3(wo, h);
    V3 wi = reflect3(wo, h);
    if (wi.z <= 0.0f || vh <= 0.0f) return wi;
    *pdf = ggx_ndf(h.z, rough) * h.z / (4.0f * vh);
    return wi;
}

/* Estimator weight w = f(wi,wo) cos_i / p(wi).
   Substituting the Cook-Torrance BRDF and p(wi) = D nh / (4 vh) gives
       w = F G (v.h) / (cos_o (n.h)).
   rough = 0 limit (G -> 1, h -> n): perfect mirror, weight = F. */
static float ggx_estimator_weight(V3 wi, V3 wo, float rough, float f0) {
    float ci = wi.z, co = wo.z;
    if (ci <= 0.0f || co <= 0.0f) return 0.0f;
    if (rough <= 0.0f) {
        V3 n = { 0.0f, 0.0f, 1.0f };
        V3 mirror = reflect3(wo, n);
        if (fabsf(wi.x - mirror.x) > 1e-4f ||
            fabsf(wi.y - mirror.y) > 1e-4f ||
            fabsf(wi.z - mirror.z) > 1e-4f) return 0.0f;
        return fresnel_schlick(co, f0);
    }
    V3 h = norm3(add3(wi, wo));
    float nh = h.z, vh = dot3(wo, h);
    if (nh <= 0.0f || vh <= 0.0f) return 0.0f;
    float g = smith_g(ci, co, rough * rough);
    float f = fresnel_schlick(vh, f0);
    return f * g * vh / (co * nh);
}

/* deterministic LCG so the check output is reproducible */
static uint32_t rng_state = 0x12345678u;
static float rng01(void) {
    rng_state = rng_state * 1664525u + 1013904223u;
    return (float)(rng_state >> 8) / 16777216.0f;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* quadrature of integral D cos d omega = 2 pi int D(theta) cos sin d theta */
static double pdf_quadrature(float rough) {
    const int K = 20000;
    double dt = (PI / 2.0) / K, acc = 0.0;
    for (int k = 0; k < K; k++) {
        double th = (k + 0.5) * dt;
        float nh = (float)cos(th);
        acc += (double)ggx_ndf(nh, rough) * (double)nh * sin(th) * dt;
    }
    return 2.0 * PI * acc;
}

/* empirical density / CDF of sampled half vectors vs analytic p_h.
   The narrow-lobe case cannot be checked with midpoint-density bins
   (a bin wider than the lobe makes the midpoint rule meaningless), so
   the strict assertion compares each bin's count against the analytic
   CDF difference; the density table is printed for the resolved case. */
static void pdf_histogram_check(float rough, int n, int bins) {
    long cnt[64];
    memset(cnt, 0, sizeof(cnt));
    double dt = (PI / 2.0) / bins;
    for (int i = 0; i < n; i++) {
        V3 h = ggx_sample_half(rng01(), rng01(), rough);
        double th = acos(fmin(fmax((double)h.z, 0.0), 1.0));
        int b = (int)(th / dt);
        if (b >= bins) b = bins - 1;
        cnt[b]++;
    }
    double a2 = (double)rough * rough * rough * rough;  /* (rough^2)^2 */
    float peak = ggx_ndf(1.0f, rough);   /* p_h peaks at theta = 0 */
    double max_rel = 0.0, chi2 = 0.0;
    double cum_emp = 0.0, max_cdf_err = 0.0;
    printf("  theta_lo theta_hi   pdf_emp   pdf_analytic  rel_err  count/sigma\n");
    for (int b = 0; b < bins; b++) {
        double lo = b * dt, hi = (b + 1) * dt, mid = 0.5 * (lo + hi);
        double omega = 2.0 * PI * (cos(lo) - cos(hi));
        double emp = (double)cnt[b] / ((double)n * omega);
        float nh = (float)cos(mid);
        double ana = (double)ggx_ndf(nh, rough) * (double)nh;
        double rel = ana > 1e-12 ? fabs(emp - ana) / ana : 0.0;
        if (ana > 0.05 * (double)peak) max_rel = fmax(max_rel, rel);
        /* expected bin count from the analytic marginal CDF;
           Pearson chi-square stays well-behaved for sparse bins */
        double f_lo = b == 0 ? 0.0 :
                      (tan(lo) * tan(lo)) / (a2 + tan(lo) * tan(lo));
        double f_hi = (tan(hi) * tan(hi)) / (a2 + tan(hi) * tan(hi));
        double exp_cnt = (double)n * (f_hi - f_lo);
        double sigma = exp_cnt >= 1.0 ?
                       fabs((double)cnt[b] - exp_cnt) / sqrt(exp_cnt) : 0.0;
        if (exp_cnt >= 1.0) chi2 += ((double)cnt[b] - exp_cnt) *
                                    ((double)cnt[b] - exp_cnt) / exp_cnt;
        cum_emp += (double)cnt[b] / (double)n;
        max_cdf_err = fmax(max_cdf_err, fabs(cum_emp - f_hi));
        if (b % 6 == 0)
            printf("  %8.3f %8.3f %10.4f %12.4f %8.4f %10.2f\n",
                   lo * 180.0 / PI, hi * 180.0 / PI, emp, ana, rel, sigma);
    }
    printf("  pearson_chi2=%.1f / dof=%d (tol %.1f), max_cdf_err=%.5f"
           " (tol %.5f)\n", chi2, bins, 2.0 * bins,
           max_cdf_err, 4.0 / sqrt((double)n));
    check("hist_count_matches_cdf", chi2 < 2.0 * bins);
    check("hist_cdf_matches_analytic", max_cdf_err < 4.0 / sqrt((double)n));
    if (rough >= 0.30f) {
        printf("  (resolved bins) max_density_rel_err=%.4f (tol 0.08)\n",
               max_rel);
        check("hist_pdf_matches_analytic", max_rel < 0.08);
    }
}

typedef struct { double mean; double var; } MStat;

/* one estimator: method 0 = cosine sampling, 1 = GGX sampling.
   returns mean and variance of the per-trial integral estimates */
static MStat mc_integral(int method, V3 wo, float rough, float f0,
                         int spp, int trials) {
    double *est = (double*)malloc(sizeof(double) * (size_t)trials);
    for (int t = 0; t < trials; t++) {
        double sum = 0.0;
        for (int s = 0; s < spp; s++) {
            float w = 0.0f;
            if (method == 0) {
                V3 wi = cosine_sample(rng01(), rng01());
                if (wi.z > 0.0f) {
                    float f = brdf_ggx(wi, wo, rough, f0);
                    w = f * PI;                       /* f cos / (cos/pi) */
                }
            } else {
                float pdf;
                V3 wi = ggx_sample_wi(wo, rng01(), rng01(), rough, &pdf);
                if (pdf > 0.0f)
                    w = brdf_ggx(wi, wo, rough, f0) * wi.z / pdf;
            }
            sum += (double)w;
        }
        est[t] = sum / spp;
    }
    MStat m;
    double acc = 0.0;
    for (int t = 0; t < trials; t++) acc += est[t];
    m.mean = acc / trials;
    acc = 0.0;
    for (int t = 0; t < trials; t++) acc += (est[t] - m.mean) * (est[t] - m.mean);
    m.var = acc / (trials - 1);
    free(est);
    return m;
}

int main(void) {
    check("lambert_unit", fabsf(brdf_lambert(PI) - 1.0f) < 1e-4f);
    check("lambert_zero", brdf_lambert(0.0f) == 0.0f);

    check("ggx_peak", ggx_ndf(1.0f, 0.5f) > ggx_ndf(0.2f, 0.5f));
    check("ggx_positive", ggx_ndf(0.5f, 0.3f) > 0.0f);

    check("pdf_half", fabsf(cosine_pdf(0.5f) - 0.5f / PI) < 1e-6f);

    V3 s = cosine_sample(0.5f, 0.25f);
    check("sample_unit", fabsf(sqrtf(dot3(s, s)) - 1.0f) < 1e-4f);
    check("sample_z_nonneg", s.z >= 0.0f);

    V3 up = { 0, 0, 1 };
    check("sample_hemisphere", dot3(s, up) >= -1e-6f);

    /* mean of cos over samples ~ 2/3 for cosine sampling */
    double acc = 0.0; int n = 512;
    for (int i = 0; i < n; i++) {
        float u1 = (float)((i % 16) + 0.5f) / 16.0f;
        float u2 = (float)((i / 16) + 0.5f) / 32.0f;
        V3 v = cosine_sample(u1, u2);
        acc += v.z;
    }
    check("mean_cos_close", fabs(acc / n - 2.0 / 3.0) < 0.02);

    /* ---- GGX pdf normalization: integral of p_h must be 1 ---- */
    check("ggx_norm_r025", fabs(pdf_quadrature(0.25f) - 1.0) < 1e-3);
    check("ggx_norm_r050", fabs(pdf_quadrature(0.50f) - 1.0) < 1e-3);
    check("ggx_norm_r090", fabs(pdf_quadrature(0.90f) - 1.0) < 1e-3);

    /* ---- sampled half vectors vs analytic pdf, pointwise ---- */
    printf("ggx half-vector histogram (rough=0.40, n=300000):\n");
    pdf_histogram_check(0.40f, 300000, 48);
    printf("ggx half-vector CDF check (rough=0.10, n=300000):\n");
    pdf_histogram_check(0.10f, 300000, 48);

    /* ---- Jacobian: direction pdf must equal p_h / (4 wo.h) ---- */
    V3 wo45 = { 0.70710678f, 0.0f, 0.70710678f };
    int jac_ok = 1, ident_ok = 1, finite_ok = 1;
    for (int i = 0; i < 64; i++) {
        float pdf;
        V3 wi = ggx_sample_wi(wo45, rng01(), rng01(), 0.35f, &pdf);
        if (pdf <= 0.0f) continue;
        V3 h = norm3(add3(wi, wo45));
        float expect = ggx_half_pdf(h, 0.35f) / (4.0f * dot3(wo45, h));
        jac_ok &= fabsf(pdf - expect) < 2e-5f * fmaxf(1.0f, expect);
        jac_ok &= fabsf(ggx_dir_pdf(wi, wo45, 0.35f) - expect) <
                  2e-5f * fmaxf(1.0f, expect);
        /* f cos / p must equal the closed-form estimator weight */
        float w1 = brdf_ggx(wi, wo45, 0.35f, 0.04f) * wi.z / pdf;
        float w2 = ggx_estimator_weight(wi, wo45, 0.35f, 0.04f);
        ident_ok &= fabsf(w1 - w2) < 1e-3f * fmaxf(1.0f, w1);
        finite_ok &= finite_f(pdf) && finite_f(w1) && finite_f(w2);
    }
    check("ggx_jacobian_pdf", jac_ok);
    check("ggx_estimator_identity", ident_ok);
    check("ggx_all_finite", finite_ok);

    /* ---- equal-sample variance: cosine vs GGX sampling ---- */
    const int SPP = 64, TRIALS = 500;
    float roughs[3] = { 0.10f, 0.25f, 0.50f };
    printf("integral of GGX BRDF * cos over hemisphere (spp=%d, trials=%d):\n",
           SPP, TRIALS);
    printf(" rough  mean_cosine   mean_ggx   var_cosine     var_ggx"
           "     var_ratio  stddev_drop\n");
    int variance_ok = 1, mean_ok = 1;
    for (int k = 0; k < 3; k++) {
        float r = roughs[k];
        MStat c = mc_integral(0, wo45, r, 0.04f, SPP, TRIALS);
        MStat g = mc_integral(1, wo45, r, 0.04f, SPP, TRIALS);
        double ratio = c.var / g.var;
        printf(" %4.2f  %11.6f %11.6f %13.3e %13.3e %9.1fx  %7.2fx\n",
               r, c.mean, g.mean, c.var, g.var, ratio, sqrt(ratio));
        variance_ok &= ratio > 1.5;
        if (r >= 0.25f)
            mean_ok &= fabs(c.mean - g.mean) < 0.15 * fmax(c.mean, g.mean);
    }
    check("ggx_variance_reduced", variance_ok);
    check("two_methods_agree", mean_ok);

    /* ---- rough = 0: delta lobe, deterministic mirror reflection ---- */
    V3 h0 = ggx_sample_half(0.37f, 0.91f, 0.0f);
    check("zero_rough_half_is_n", h0.x == 0.0f && h0.y == 0.0f && h0.z == 1.0f);
    V3 mirror = reflect3(wo45, up);
    float pdf0;
    V3 wi0 = ggx_sample_wi(wo45, 0.37f, 0.91f, 0.0f, &pdf0);
    check("zero_rough_mirror",
          fabsf(wi0.x - mirror.x) < 1e-6f && fabsf(wi0.y) < 1e-6f &&
          fabsf(wi0.z - mirror.z) < 1e-6f);
    check("zero_rough_pdf_not_finite", pdf0 == 0.0f);
    float w0 = ggx_estimator_weight(wi0, wo45, 0.0f, 0.04f);
    check("zero_rough_weight_is_fresnel",
          finite_f(w0) && fabsf(w0 - fresnel_schlick(wo45.z, 0.04f)) < 1e-6f);
    check("zero_rough_ndf_regular", ggx_ndf(0.5f, 0.0f) == 0.0f &&
          finite_f(ggx_ndf(0.5f, 0.0f)));
    /* concentration: as rough -> 0 the sampled lobe collapses onto n */
    double mean_theta = 0.0;
    for (int i = 0; i < 10000; i++) {
        V3 h = ggx_sample_half(rng01(), rng01(), 0.02f);
        mean_theta += acos(fmin(fmax((double)h.z, 0.0), 1.0));
    }
    check("small_rough_concentrates", mean_theta / 10000.0 < 0.05);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
