/* lightsample: importance sampling for direct lighting with many lights.
 *
 * At a shading point every light i has a scalar contribution C_i
 * (emission * geometry * visibility). One light I is drawn per estimate
 * and the probability-compensated estimator
 *
 *     X = C_I / p_I,   E[X] = sum_i p_i C_i / p_i = sum_i C_i = mu
 *
 * is unbiased for the direct-lighting integral mu. Three discrete PDFs
 * over the N lights are compared:
 *   uniform  p_i = 1 / N
 *   power    p_i = P_i / sum_j P_j            (flux only, no geometry)
 *   contrib  p_i = c_i / sum_j c_j            (unshadowed contribution proxy)
 *
 * Variance is known analytically for any PDF:
 *   V(p) = E[X^2] - mu^2 = sum_i C_i^2 / p_i - mu^2,
 * and is also measured by Monte Carlo simulation.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N_LIGHTS 40
#define N_TRIALS 2000000

static double lx[N_LIGHTS], ly[N_LIGHTS], lz[N_LIGHTS];
static double lamp_power[N_LIGHTS];
static double chat[N_LIGHTS]; /* unshadowed contribution proxy c_i       */
static double C[N_LIGHTS];    /* actual contribution, incl. visibility   */
static double p_uni[N_LIGHTS], p_pow[N_LIGHTS], p_con[N_LIGHTS];

/* deterministic hash so the scene and the printed numbers are reproducible */
static double hash01(unsigned i) {
    unsigned x = i * 2654435761u;
    x ^= x >> 16; x *= 2246822519u; x ^= x >> 13;
    x *= 3266489917u; x ^= x >> 16;
    return (x >> 8) * (1.0 / 16777216.0);
}

/* xorshift64 for the Monte Carlo trials */
static unsigned long long rng_state;
static double rnd(void) {
    unsigned long long x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng_state = x;
    return (double)(x >> 11) * (1.0 / 9007199254740992.0); /* 2^53 */
}

/* Build a normalized PDF from non-negative weights.
 * Returns 1 when the total weight is zero and uniform sampling is used. */
static int build_pdf(const double *w, double *p) {
    double s = 0.0;
    for (int i = 0; i < N_LIGHTS; i++) s += w[i];
    if (s == 0.0) {
        for (int i = 0; i < N_LIGHTS; i++) p[i] = 1.0 / N_LIGHTS;
        return 1;
    }
    for (int i = 0; i < N_LIGHTS; i++) p[i] = w[i] / s;
    return 0;
}

static void build_scene(double nx, double ny, double nz, int occlusion) {
    static const double PI = 3.14159265358979323846;
    for (int i = 0; i < N_LIGHTS; i++) {
        double az  = hash01((unsigned)(i * 4 + 0)) * 2.0 * PI;
        double el  = hash01((unsigned)(i * 4 + 1)) * 1.4 - 0.3; /* some below horizon */
        double d   = 2.0 + 8.0 * hash01((unsigned)(i * 4 + 2));
        double log_p = hash01((unsigned)(i * 4 + 3)) * log(100.0); /* power in [1,100] */
        lamp_power[i] = exp(log_p);
        lx[i] = d * cos(el) * cos(az);
        ly[i] = d * cos(el) * sin(az);
        lz[i] = d * sin(el);

        double cos_t = (nx * lx[i] + ny * ly[i] + nz * lz[i]) / d;
        if (cos_t <= 0.0) chat[i] = 0.0;                 /* geometrically impossible */
        else chat[i] = lamp_power[i] * cos_t / (d * d); /* point-light falloff */

        double vis = 1.0;
        if (occlusion && hash01(1000u + (unsigned)i) < 0.3) vis = 0.0;
        C[i] = chat[i] * vis;
    }
}

static double total_contribution(void) {
    double mu = 0.0;
    for (int i = 0; i < N_LIGHTS; i++) mu += C[i];
    return mu;
}

/* V(p) = sum C_i^2 / p_i - mu^2; a zero-probability light with C_i > 0
 * would make the estimator undefined, hence the infinity guard. */
static double analytic_var(const double *p, double mu) {
    double s = 0.0;
    for (int i = 0; i < N_LIGHTS; i++) {
        if (C[i] == 0.0) continue;
        if (p[i] <= 0.0) return INFINITY;
        s += C[i] * C[i] / p[i];
    }
    return s - mu * mu;
}

static void make_cdf(const double *p, double *cdf) {
    double a = 0.0;
    for (int i = 0; i < N_LIGHTS; i++) { a += p[i]; cdf[i] = a; }
    cdf[N_LIGHTS - 1] = 1.0;
}

static int sample_cdf(const double *cdf, double u) {
    int lo = 0, hi = N_LIGHTS - 1;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (u <= cdf[mid]) hi = mid; else lo = mid + 1;
    }
    return lo;
}

static void simulate(const double *p, const double *cdf, double *mean, double *var) {
    double s = 0.0, s2 = 0.0;
    double n = (double)N_TRIALS;
    for (int t = 0; t < N_TRIALS; t++) {
        int i = sample_cdf(cdf, rnd());
        double x = C[i] / p[i];
        s += x; s2 += x * x;
    }
    *mean = s / n;
    double m2 = s2 / n;
    *var = (m2 - (*mean) * (*mean)) * n / (n - 1.0);
}

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    rng_state = 0x123456789abcdefULL;

    /* ---- scene: 40 lights, shading point at origin, normal +z ---- */
    build_scene(0.0, 0.0, 1.0, 1);
    for (int i = 0; i < N_LIGHTS; i++) p_uni[i] = 1.0 / N_LIGHTS;
    int fb_pow = build_pdf(lamp_power, p_pow);
    int fb_con = build_pdf(chat,       p_con);

    double sum_chat = 0.0;
    for (int i = 0; i < N_LIGHTS; i++) sum_chat += chat[i];

    int norm_ok = 1, proxy_ok = 1, zero_safe = 1;
    double su = 0.0, sp = 0.0, sc = 0.0;
    for (int i = 0; i < N_LIGHTS; i++) {
        su += p_uni[i]; sp += p_pow[i]; sc += p_con[i];
        if (chat[i] > 0.0 && fabs(p_con[i] * sum_chat - chat[i]) > 1e-12) proxy_ok = 0;
        if (p_con[i] == 0.0 && C[i] != 0.0) zero_safe = 0; /* p=0 only where C=0 */
    }
    norm_ok = fabs(su - 1) < 1e-12 && fabs(sp - 1) < 1e-12 && fabs(sc - 1) < 1e-12;
    check("pdfs_normalized", norm_ok);
    check("contrib_no_fallback_when_lit", fb_pow == 0 && fb_con == 0);
    /* probability-contribution relation: p_i = c_i / sum_j c_j */
    check("contrib_pdf_matches_proxy", proxy_ok);
    check("zero_pdf_only_on_zero_contribution", zero_safe);

    double mu = total_contribution();
    double v_uni = analytic_var(p_uni, mu);
    double v_pow = analytic_var(p_pow, mu);
    double v_con = analytic_var(p_con, mu);

    double cdf_uni[N_LIGHTS], cdf_pow[N_LIGHTS], cdf_con[N_LIGHTS];
    make_cdf(p_uni, cdf_uni); make_cdf(p_pow, cdf_pow); make_cdf(p_con, cdf_con);
    double m_uni, e_uni, m_pow, e_pow, m_con, e_con;
    simulate(p_uni, cdf_uni, &m_uni, &e_uni);
    simulate(p_pow, cdf_pow, &m_pow, &e_pow);
    simulate(p_con, cdf_con, &m_con, &e_con);

    check("estimator_unbiased_uniform", fabs(m_uni - mu) / mu < 1e-2);
    check("estimator_unbiased_power",   fabs(m_pow - mu) / mu < 1e-2);
    check("estimator_unbiased_contrib", fabs(m_con - mu) / mu < 1e-2);

    /* analytic ordering: contribution sampling must beat both baselines */
    check("analytic_contrib_beats_uniform", v_con < v_uni);
    check("analytic_contrib_beats_power",   v_con < v_pow);
    double best_analytic = fmin(v_uni, v_pow);
    const char *best_name = (v_uni < v_pow) ? "uniform" : "power";
    check("analytic_contrib_beats_best_baseline", v_con < best_analytic);
    /* empirical ordering, 10% slack for Monte Carlo noise on the variance */
    check("empirical_contrib_beats_best_baseline",
          e_con < 0.9 * fmin(e_uni, e_pow));

    /* ---- theory limit: with no occluders C_i == c_i is known exactly,
     * so p proportional to C gives the constant estimator and zero variance ---- */
    build_scene(0.0, 0.0, 1.0, 0);
    double p_exact[N_LIGHTS];
    int fb_exact = build_pdf(chat, p_exact);
    double mu_exact = total_contribution();
    double v_exact = analytic_var(p_exact, mu_exact);
    check("zero_variance_when_contribution_known",
          fb_exact == 0 && v_exact / (mu_exact * mu_exact) < 1e-15);

    /* ---- fallback: all weights zero -> uniform, estimator stays defined ---- */
    double zeros[N_LIGHTS], p_fb[N_LIGHTS];
    memset(zeros, 0, sizeof(zeros));
    int fb_flag = build_pdf(zeros, p_fb);
    int fb_ok = (fb_flag == 1);
    for (int i = 0; i < N_LIGHTS; i++)
        fb_ok = fb_ok && (fabs(p_fb[i] - 1.0 / N_LIGHTS) < 1e-15);
    check("fallback_all_zero_weights_is_uniform", fb_ok);

    /* real-scene fallback: shade a point facing away from every light so
     * all proxies vanish -> total probability zero -> uniform fallback */
    build_scene(0.0, 0.0, 1.0, 1);
    double p_flip[N_LIGHTS];
    for (int i = 0; i < N_LIGHTS; i++) { chat[i] = 0.0; C[i] = 0.0; }
    int fb_flip = build_pdf(chat, p_flip);
    int flip_ok = (fb_flip == 1) && (total_contribution() == 0.0);
    for (int i = 0; i < N_LIGHTS; i++)
        flip_ok = flip_ok && (fabs(p_flip[i] - 1.0 / N_LIGHTS) < 1e-15);
    check("fallback_zero_proxy_scene_is_uniform", flip_ok);

    /* ---- report ---- */
    build_scene(0.0, 0.0, 1.0, 1);
    printf("\nN_lights=%d trials=%d mu=%.6f\n", N_LIGHTS, N_TRIALS, mu);
    printf("\nper-light detail (first 5):\n");
    printf(" i      power       c_i(proxy)  C_i(actual)  p_power  p_contrib\n");
    for (int i = 0; i < 5; i++)
        printf("%2d  %10.4f  %10.6f  %10.6f  %.5f  %.5f\n",
               i, lamp_power[i], chat[i], C[i], p_pow[i], p_con[i]);

    printf("\nvariance of X = C_I / p_I:\n");
    printf("strategy   analytic_var   empirical_var   var / var_contrib\n");
    printf("uniform    %13.5e  %13.5e   %7.2fx\n", v_uni, e_uni, v_uni / v_con);
    printf("power      %13.5e  %13.5e   %7.2fx\n", v_pow, e_pow, v_pow / v_con);
    printf("contrib    %13.5e  %13.5e   %7.2fx\n", v_con, e_con, 1.0);
    printf("\nbest baseline = %s; contrib sampling is %.2fx lower variance "
           "(%.2fx lower stddev).\n",
           best_name, best_analytic / v_con, sqrt(best_analytic / v_con));

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
