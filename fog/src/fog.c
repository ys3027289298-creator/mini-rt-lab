/* fog: free-path sampling and transmittance estimation in participating
   media.

   Scattering distances follow the exponential free-path distribution.  In a
   homogeneous medium the density compensation is the analytic attenuation
   weight exp(-sigma*t); in a heterogeneous medium weighted delta tracking
   uses a majorant and compensates every tentative collision by
   (1 - sigma(t)/sigma_bar).  Tests compare sampled attenuation against the
   analytic transmittance with explicit Monte Carlo error bounds, and check
   that noise shrinks as the sample count grows.  Zero, negative, or non-finite
   density is rejected instead of producing an infinite or negative distance. */
#include <stdio.h>
#include <math.h>
#include <string.h>

/* deterministic xorshift64* so every run reports the same numbers */
typedef struct { unsigned long long state; } Rng;

static unsigned long long rng_u64(Rng *r) {
    unsigned long long x = r->state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    r->state = x;
    return x * 2685821657736338717ULL;
}

/* uniform in [0, 1) */
static double rng_real(Rng *r) {
    return (double)(rng_u64(r) >> 11) * (1.0 / 9007199254740992.0);
}

/* Free-path sampling: t ~ sigma * exp(-sigma*t), t >= 0.
   The density must be finite and strictly positive; zero or negative density
   has no valid exponential distribution and is rejected. */
static int sample_free_path(double sigma, double u, double *t) {
    if (!isfinite(sigma) || sigma <= 0.0) return 0;
    if (!(u >= 0.0 && u < 1.0)) return 0;
    *t = -log1p(-u) / sigma;
    return 1;
}

/* Homogeneous medium: sample a distance and compensate its weight by the
   density along the path, exp(-sigma*t).  E[weight] = 1/2 analytically. */
static int sample_attenuation(double sigma, double u, double *t, double *weight) {
    if (!sample_free_path(sigma, u, t)) return 0;
    *weight = exp(-sigma * *t);
    return 1;
}

/* Heterogeneous fog with exponential density profile
   sigma(t) = sigma0 * exp(-t / height), sigma0 > 0, height > 0. */
typedef struct { double sigma0; double height; } Fog;

static int fog_valid(const Fog *f) {
    return f && isfinite(f->sigma0) && f->sigma0 > 0.0 &&
           isfinite(f->height) && f->height > 0.0;
}

static int fog_density(const Fog *f, double t, double *sigma) {
    if (!fog_valid(f) || !isfinite(t) || t < 0.0) return 0;
    *sigma = f->sigma0 * exp(-t / f->height);
    return 1;
}

/* Analytic optical depth tau(t) = int_0^t sigma(s) ds. */
static double fog_optical_depth(const Fog *f, double t) {
    return f->sigma0 * f->height * (1.0 - exp(-t / f->height));
}

/* Weighted delta tracking: unbiased estimate of exp(-tau(tmax)) using
   tentative collisions sampled at the majorant rate sigma_bar >= max sigma.
   At each tentative collision the weight is multiplied by
   1 - sigma(t)/sigma_bar, which compensates for the local density. */
static int fog_transmittance_delta(const Fog *f, double sigma_bar, double tmax,
                                   Rng *rng, double *transmittance) {
    if (!fog_valid(f)) return 0;
    if (!isfinite(sigma_bar) || sigma_bar <= 0.0 || sigma_bar < f->sigma0) return 0;
    if (!isfinite(tmax) || tmax < 0.0) return 0;

    double weight = 1.0, t = 0.0;
    int done = 0;
    for (int i = 0; i < 1000000; i++) {
        t += -log1p(-rng_real(rng)) / sigma_bar;
        if (t >= tmax) { done = 1; break; }
        double sigma;
        if (!fog_density(f, t, &sigma)) return 0;
        if (sigma > sigma_bar) return 0;   /* majorant violated */
        weight *= 1.0 - sigma / sigma_bar;
        if (weight <= 0.0) { *transmittance = 0.0; return 1; }
    }
    if (!done) return 0;
    *transmittance = weight;
    return 1;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* finite-sample bound: 5 standard errors of a binomial proportion */
static double binomial_bound(double p, int n) {
    double se = sqrt(p * (1.0 - p) / (double)n);
    return 5.0 * se;
}

int main(void) {
    double t = 0.0, w = 0.0;

    /* zero / negative / non-finite density is rejected */
    check("reject_zero_density", !sample_free_path(0.0, 0.5, &t));
    check("reject_negative_density", !sample_free_path(-1.0, 0.5, &t));
    check("reject_nan_density", !sample_free_path(NAN, 0.5, &t));
    check("reject_inf_density", !sample_free_path(INFINITY, 0.5, &t));

    Fog zero_fog = { 0.0, 2.0 };
    Fog neg_fog = { -1.0, 2.0 };
    Fog neg_height = { 1.0, -2.0 };
    check("reject_zero_profile", !fog_valid(&zero_fog));
    check("reject_negative_profile", !fog_valid(&neg_fog));
    check("reject_negative_height", !fog_valid(&neg_height));

    double dummy;
    Rng bad_rng = { 88172645463325252ULL };
    check("reject_zero_delta_density",
          !fog_transmittance_delta(&zero_fog, 1.0, 1.0, &bad_rng, &dummy));
    check("reject_negative_delta_density",
          !fog_transmittance_delta(&neg_fog, 1.0, 1.0, &bad_rng, &dummy));

    /* homogeneous free paths: mean and CDF against the analytic exponential */
    const double sigma = 1.7;
    const int n = 50000;
    const int nx = 5;
    const double xs[5] = { 0.1, 0.25, 0.5, 1.0, 2.0 };
    double depth = 0.9;
    int counts[5] = { 0 };
    int survived = 0;
    double sum_t = 0.0, sum_w = 0.0, sum_w2 = 0.0;
    int samples_ok = 1;

    Rng rng = { 88172645463325252ULL };
    for (int i = 0; i < n; i++) {
        double u = rng_real(&rng);
        if (!sample_attenuation(sigma, u, &t, &w) || !(t >= 0.0) || !(w > 0.0 && w <= 1.0)) {
            samples_ok = 0;
            break;
        }
        sum_t += t;
        sum_w += w;
        sum_w2 += w * w;
        for (int j = 0; j < nx; j++) if (t <= xs[j]) counts[j]++;
        if (t > depth) survived++;
    }

    check("samples_in_range", samples_ok);
    check("free_path_mean", samples_ok && fabs(sum_t / n - 1.0 / sigma) < 0.03 / sigma);

    /* E[exp(-sigma*t)] for t ~ Exp(sigma) is sigma/(2*sigma) = 1/2 */
    check("attenuation_weight_mean", samples_ok && fabs(sum_w / n - 0.5) < 0.01);

    /* empirical CDF vs 1 - exp(-sigma*x), each within 5 standard errors */
    int cdf_ok = samples_ok;
    for (int j = 0; j < nx && cdf_ok; j++) {
        double f = 1.0 - exp(-sigma * xs[j]);
        double emp = (double)counts[j] / (double)n;
        if (fabs(emp - f) > binomial_bound(f, n)) cdf_ok = 0;
    }
    check("free_path_cdf_bound", cdf_ok);

    /* sampled survival probability vs analytic transmittance at fixed depth */
    double trans = exp(-sigma * depth);
    double emp_trans = (double)survived / (double)n;
    check("homogeneous_transmittance_bound",
          samples_ok && fabs(emp_trans - trans) < binomial_bound(trans, n));

    /* heterogeneous fog: weighted delta tracking vs analytic transmittance */
    Fog fog = { 1.0, 2.0 };
    double tmax = 4.0;
    double T_an = exp(-fog_optical_depth(&fog, tmax));
    const int nd = 200000;
    double sum_T = 0.0, sum_T2 = 0.0, est = 0.0;
    int delta_ok = 1;

    Rng drng = { 1442695040888963407ULL };
    for (int i = 0; i < nd; i++) {
        double T;
        if (!fog_transmittance_delta(&fog, fog.sigma0, tmax, &drng, &T) ||
            !(T >= 0.0 && T <= 1.0)) {
            delta_ok = 0;
            break;
        }
        sum_T += T;
        sum_T2 += T * T;
    }
    if (delta_ok) est = sum_T / nd;
    double est_var = delta_ok ? sum_T2 / nd - est * est : 0.0;
    double est_se = delta_ok && est_var > 0.0 ? sqrt(est_var / nd) : 0.0;
    double rel_err = delta_ok ? fabs(est - T_an) / T_an : 1.0;

    check("delta_weight_in_unit", delta_ok);
    check("delta_transmittance_bound", delta_ok && rel_err < 0.01);
    check("delta_transmittance_se_bound", delta_ok && fabs(est - T_an) <= 5.0 * est_se);

    /* more samples -> less noise: RMSE over independent batches falls */
    const int batches = 64, n_small = 500, n_large = 2000;
    double sse_small = 0.0, sse_large = 0.0;
    int conv_ok = 1;
    for (int b = 0; b < batches && conv_ok; b++) {
        Rng r1 = { 0x9e3779b97f4a7c15ULL + (unsigned long long)b * 2ULL };
        Rng r2 = { 0x9e3779b97f4a7c15ULL + (unsigned long long)b * 2ULL + 1ULL };
        double s1 = 0.0, s2 = 0.0;
        for (int i = 0; i < n_large; i++) {
            double T;
            if (!fog_transmittance_delta(&fog, fog.sigma0, tmax, &r1, &T)) { conv_ok = 0; break; }
            if (i < n_small) s1 += T;
            if (!fog_transmittance_delta(&fog, fog.sigma0, tmax, &r2, &T)) { conv_ok = 0; break; }
            s2 += T;
        }
        if (!conv_ok) break;
        double e1 = s1 / n_small - T_an;
        double e2 = s2 / n_large - T_an;
        sse_small += e1 * e1;
        sse_large += e2 * e2;
    }
    double rmse_small = sqrt(sse_small / batches);
    double rmse_large = sqrt(sse_large / batches);
    /* 4x samples should cut RMSE about in half; allow headroom for batch noise */
    check("noise_converges", conv_ok && rmse_large < 0.75 * rmse_small);

    printf("delta_estimate=%.6f delta_analytic=%.6f rmse_500=%.6f rmse_2000=%.6f\n",
           est, T_an, rmse_small, rmse_large);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
