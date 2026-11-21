
/* sample: stratified and low-discrepancy sampling with variance estimate,
   plus error confidence for iid Monte Carlo estimators:
     - confidence interval for the estimate, from the sample variance
     - sample count needed to reach a target standard error
     - measured-error / theoretical-error ratio, checked on a convergence
       curve, with an upper bound on the sample count implied by the
       uncertainty of the variance estimate itself.
   The iid assumption matters: the stratified estimate below reuses each
   strata once, so its samples are not independent and the usual
   sigma/sqrt(n) formula does not apply to it.  Everything here is for
   plain iid sampling unless stated otherwise. */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define N_MAX 256
#define MC_MAX 4096

typedef struct { double mean; double var; int n; } Stat;

/* Radical inverse of i in the given base, as a number in [0,1). */
static double halton2(unsigned i, unsigned base) {
    double f = 1.0, r = 0.0; unsigned k = i;
    while (k > 0) { f /= base; r += f * (k % base); k /= base; }
    return r;
}

/* estimate integral of f over [0,1) with stratified samples.
   Each stratum is used exactly once, so these samples are *not*
   independent: the iid machinery below does not apply to this mean. */
static Stat stratify(int strata, double (*f)(double)) {
    Stat s; memset(&s, 0, sizeof(s));
    double sum = 0.0, sum2 = 0.0;
    int n = 0;
    for (int i = 0; i < strata; i++) {
        double u = ((double)i + 0.5) / (double)strata;
        double v = f(u);
        sum += v; sum2 += v * v; n++;
    }
    s.n = n;
    s.mean = sum / n;
    double var = (n > 1) ? (sum2 / n - s.mean * s.mean) * n / (n - 1) : 0.0;
    s.var = var > 0 ? var : 0;
    return s;
}

/* ------------------------------------------------------------------ */
/* Error confidence for iid samples                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    int n;
    double mean, var;   /* var: unbiased sample variance, sum (x-mean)^2/(n-1) */
    double se;          /* sqrt(var / n), the estimated standard error         */
    double lo, hi;      /* mean +- z * se                                       */
    int valid;          /* 0 when no interval can be formed                     */
} McStat;

static double mc_mean(const double* x, int n) {
    if (n <= 0) return NAN;
    double s = 0.0;
    for (int i = 0; i < n; i++) s += x[i];
    return s / (double)n;
}

static double mc_var(const double* x, int n) {
    if (n < 2) return 0.0;
    double m = mc_mean(x, n), s = 0.0;
    for (int i = 0; i < n; i++) { double d = x[i] - m; s += d * d; }
    return s / (double)(n - 1);
}

/* Standard error of the mean.  Zero variance or fewer than two samples
   gives zero: no spread was observed, so there is nothing to quote. */
static double mc_std_error(double var, int n) {
    if (n < 2 || !(var > 0.0)) return 0.0;
    return sqrt(var / (double)n);
}

/* Inverse standard normal CDF (Acklam's rational approximation, relative
   error about 1e-9).  We use the normal quantile rather than Student-t:
   no t quantile exists in libm, and t -> z as n grows.  At small n the
   interval below is therefore a bit too narrow, and since sigma is
   estimated from the same samples the interval is only as good as that
   estimate.  Both caveats are checked by the coverage test in main. */
static double z_quantile(double p) {
    static const double a[6] = { -3.969683028665376e+01, 2.209460984245205e+02,
                                 -2.759285104469687e+02, 1.383577518672690e+02,
                                 -3.066479806614716e+01, 2.506628277459239e+00 };
    static const double b[5] = { -5.447609879822406e+01, 1.615858368580409e+02,
                                 -1.556989798598866e+02, 6.680131188771972e+01,
                                 -1.328068155288572e+01 };
    static const double c[6] = { -7.784894002430293e-03, -3.223964580411365e-01,
                                 -2.400758277161838e+00, -2.549732539343734e+00,
                                 4.374664141464968e+00, 2.938163982698783e+00 };
    static const double d[4] = { 7.784695709041462e-03, 3.224671290700398e-01,
                                 2.445134137142996e+00, 3.754408661907416e+00 };
    double q, r;
    if (!(p > 0.0 && p < 1.0)) return NAN;
    if (p > 0.5) return -z_quantile(1.0 - p);
    if (p < 0.02425) {
        q = sqrt(-2.0 * log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
    q = p - 0.5; r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
}

/* Two-sided normal critical value for coverage 1 - alpha. */
static double z_crit(double alpha) {
    if (!(alpha > 0.0 && alpha < 1.0)) return NAN;
    return z_quantile(1.0 - alpha / 2.0);
}

/* Confidence interval for the mean.  Returns 0 (and leaves lo/hi NaN) if
   alpha is outside (0,1) or there are fewer than two samples; a zero
   variance is accepted and yields the point interval [mean, mean]. */
static int mc_ci(const double* x, int n, double alpha, McStat* out) {
    memset(out, 0, sizeof(*out));
    out->n = n;
    out->mean = (n > 0) ? mc_mean(x, n) : NAN;
    out->var = mc_var(x, n);
    out->se = mc_std_error(out->var, n);
    out->lo = out->hi = NAN;
    double z = z_crit(alpha);
    if (isnan(z) || n < 2) return 0;
    out->lo = out->mean - z * out->se;
    out->hi = out->mean + z * out->se;
    out->valid = 1;
    return 1;
}

/* Sample count for standard error <= target, given the population
   variance.  Returns -1 when the target is unreachable: a non-positive
   variance says the estimator cannot get closer, a non-positive target
   asks for exactly zero error, and non-finite inputs have no answer.
   Never guess a finite number for those cases. */
static int mc_n_needed(double var, double target) {
    if (!(var > 0.0) || !(target > 0.0)) return -1;
    if (!isfinite(var) || !isfinite(target)) return -1;
    double n = var / (target * target);
    if (!isfinite(n) || n >= 2147483647.0) return -1;
    if (n < 2.0) n = 2.0;
    return (int)ceil(n);
}

/* One-sided upper bound on the population variance at confidence
   1 - alpha, from the uncertainty of the variance estimate:
   Var(s^2) = (mu4 - (n-3)/(n-1) * sigma^4) / n, with mu4 and sigma^4
   estimated from the same samples.  Asymptotic (normal-ish) in n; it is
   what turns a point estimate of n* into a "no more than this" count. */
static double mc_var_upper(const double* x, int n, double alpha) {
    if (n < 4) return NAN;
    double m = mc_mean(x, n), var = mc_var(x, n);
    double s4 = 0.0;
    for (int i = 0; i < n; i++) { double d = x[i] - m, d2 = d * d; s4 += d2 * d2; }
    double m4 = s4 / (double)n;
    double vs = (m4 - ((double)(n - 3) / (double)(n - 1)) * var * var) / (double)n;
    if (!(vs > 0.0)) return var;          /* no usable fourth moment: no inflation */
    double z = z_quantile(1.0 - alpha);   /* one-sided */
    if (isnan(z)) return NAN;
    return var + z * sqrt(vs);
}

/* ------------------------------------------------------------------ */
/* Self tests                                                          */
/* ------------------------------------------------------------------ */

static double f_quad(double u) { return u * u; }
static double f_const(double u) { (void)u; return 1.0; }

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Deterministic RNG so the statistical tests below are reproducible. */
static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static uint64_t splitmix64(void) {
    uint64_t z = (rng_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double urand(void) {
    return (double)(splitmix64() >> 11) * (1.0 / 9007199254740992.0);
}

#define REP_LIMIT 400
#define GRID_LEN 5

int main(void) {
    check("halton_in_unit", halton2(1, 2) >= 0.0 && halton2(1, 2) < 1.0);
    check("halton_first", fabs(halton2(0, 2)) < 1e-12);

    Stat a = stratify(64, f_quad);
    check("quad_mean_close", fabs(a.mean - 1.0 / 3.0) < 1e-3);

    Stat b = stratify(16, f_quad);
    Stat c = stratify(64, f_quad);
    check("variance_drops", c.var <= b.var + 1e-6);

    Stat d = stratify(32, f_const);
    check("const_exact", fabs(d.mean - 1.0) < 1e-12 && d.var < 1e-12);

    check("halton_distinct", halton2(1, 2) != halton2(2, 2));
    check("halton_monotone_base", halton2(1, 2) > halton2(0, 2));

    /* ---- quantiles ---- */
    check("z_crit_95", fabs(z_crit(0.05) - 1.959964) < 1e-5);
    check("z_crit_symmetric", fabs(z_crit(0.10) - 1.644854) < 1e-5);
    check("z_crit_bad_alpha", isnan(z_crit(0.0)) && isnan(z_crit(1.0)) &&
                              isnan(z_crit(-0.1)) && isnan(z_crit(NAN)));

    /* ---- zero variance: interval is a point, sample count unreachable ---- */
    double cst[8];
    for (int i = 0; i < 8; i++) cst[i] = 2.0;
    McStat cs;
    check("ci_zero_var_valid", mc_ci(cst, 8, 0.05, &cs) == 1);
    check("ci_zero_var_point", cs.se == 0.0 && cs.lo == 2.0 && cs.hi == 2.0);
    check("n_needed_zero_var_unreachable", mc_n_needed(0.0, 0.01) == -1);
    check("n_needed_zero_target_unreachable", mc_n_needed(4.0 / 45.0, 0.0) == -1);
    check("n_needed_negative_target_unreachable", mc_n_needed(4.0 / 45.0, -1.0) == -1);
    check("n_needed_nan_unreachable", mc_n_needed(NAN, 0.01) == -1 &&
                                       mc_n_needed(1.0, NAN) == -1);
    /* Fewer than two samples: no interval, and the count is still defined
       from the supplied variance (it is a planning number). */
    check("ci_one_sample_invalid", mc_ci(cst, 1, 0.05, &cs) == 0 && isnan(cs.lo));

    /* ---- theory: f(u)=u^2 on U[0,1) has mean 1/3, variance 4/45 ---- */
    check("n_needed_theory", mc_n_needed(4.0 / 45.0, 0.01) == 889);

    const int grid[GRID_LEN] = { 16, 64, 256, 1024, 4096 };
    const double true_mean = 1.0 / 3.0, true_var = 4.0 / 45.0;
    const double target = 0.01;
    const double n_star = true_var / (target * target);   /* 888.9 samples */

    static double x[MC_MAX];
    double rms[GRID_LEN] = { 0 };
    int cover = 0, upper_covers = 0, upper_valid = 0;
    double n_guess_sum = 0.0;
    int n_guess_cnt = 0;

    for (int r = 0; r < REP_LIMIT; r++) {
        for (int i = 0; i < MC_MAX; i++) {
            double u = urand();
            x[i] = u * u;
        }
        for (int g = 0; g < GRID_LEN; g++) {
            McStat s;
            mc_ci(x, grid[g], 0.05, &s);
            double err = s.mean - true_mean;
            rms[g] += err * err;
            if (grid[g] == 256) {
                if (s.lo <= true_mean && true_mean <= s.hi) cover++;
                double vu = mc_var_upper(x, grid[g], 0.05);   /* 95% one-sided */
                int nu = mc_n_needed(vu, target);
                if (nu >= 0) {
                    upper_valid++;
                    if ((double)nu >= n_star) upper_covers++;
                }
            }
            if (grid[g] == 1024) {
                int np = mc_n_needed(s.var, target);
                if (np > 0) { n_guess_sum += (double)np; n_guess_cnt++; }
            }
        }
    }

    for (int g = 0; g < GRID_LEN; g++) rms[g] = sqrt(rms[g] / (double)REP_LIMIT);

    /* measured / theoretical error.  For iid samples E[(mean-mu)^2] is
       exactly var/n, so the ratio should sit at 1 at every n, and its
       scatter over the replicates is about 1/sqrt(2*REP_LIMIT) ~ 0.035. */
    double ratio[GRID_LEN];
    int ratio_ok = 1;
    for (int g = 0; g < GRID_LEN; g++) {
        double theory = sqrt(true_var / (double)grid[g]);
        ratio[g] = rms[g] / theory;
        printf("info_n=%d info_rms=%.6f info_theory=%.6f info_ratio=%.4f\n",
               grid[g], rms[g], theory, ratio[g]);
        if (fabs(ratio[g] - 1.0) > 0.12) ratio_ok = 0;
    }
    check("measured_vs_theory_ratio", ratio_ok);

    /* interval coverage at n = 256, nominal 95% */
    printf("info_coverage=%.4f\n", (double)cover / (double)REP_LIMIT);
    check("ci_coverage_nominal", cover >= 330 && cover <= 392);   /* ~[0.825, 0.98] */

    /* convergence curve: log-log slope of rms vs n should be -1/2 */
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int g = 0; g < GRID_LEN; g++) {
        double lx = log((double)grid[g]), ly = log(rms[g]);
        sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
    }
    double slope = ((double)GRID_LEN * sxy - sx * sy) /
                   ((double)GRID_LEN * sxx - sx * sx);
    printf("info_slope=%.4f\n", slope);
    check("curve_slope_half", fabs(slope + 0.5) < 0.05);

    /* empirical sample count at the target, log-log interpolated between
       grid points, versus the theoretical 888.9. */
    double n_emp = -1.0;
    for (int g = 1; g < GRID_LEN; g++) {
        if (rms[g - 1] > target && rms[g] <= target) {
            double lx0 = log((double)grid[g - 1]), lx1 = log((double)grid[g]);
            double ly0 = log(rms[g - 1]), ly1 = log(rms[g]);
            n_emp = exp(lx0 + (ly0 - log(target)) / (ly0 - ly1) * (lx1 - lx0));
            break;
        }
    }
    double dev = (n_emp - n_star) / n_star;
    printf("info_n_star=%.1f info_n_emp=%.1f info_dev=%.3f\n", n_star, n_emp, dev);
    check("empirical_n_close", n_emp > 0.0 && fabs(dev) < 0.25);

    /* upper bound from the variance uncertainty: it must cover the true
       count for at least ~95% of replicates (one-sided 95% bound). */
    double upper_rate = upper_valid ? (double)upper_covers / (double)upper_valid : 0.0;
    printf("info_upper_cover=%.4f\n", upper_rate);
    check("n_upper_bound_covers", upper_valid == REP_LIMIT && upper_rate >= 0.92);

    /* the per-replicate planning count from the sample variance brackets
       the theoretical count on average. */
    double n_guess = n_guess_cnt ? n_guess_sum / (double)n_guess_cnt : -1.0;
    printf("info_n_from_sample_var=%.1f\n", n_guess);
    check("n_from_sample_var_close", n_guess > 0.0 && fabs(n_guess / n_star - 1.0) < 0.12);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
