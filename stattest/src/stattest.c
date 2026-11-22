
/* stattest: chi-square goodness-of-fit test with power analysis.

   Question: is the noise really random, or is the sampler biased?

   Design:
     H0: samples iid ~ p (known bin probabilities)
     H1: samples iid ~ q != p, deviation measured by chi-square divergence
         delta^2 = sum_i (q_i - p_i)^2 / p_i
     statistic  X^2 = sum_i (O_i - N p_i)^2 / (N p_i)
       under H0: X^2 ~ chi2_{K-1};  reject when X^2 > crit,
                 crit = chi2_{K-1, 1-alpha}  (false positive rate = alpha)

   Power (two-moment exact approximation):
     First-order theory says X^2 under H1 is noncentral chi2 with
     lambda = N delta^2.  That approximation OVERESTIMATES power at
     moderate lambda (measured here: 0.945 vs 0.903 empirical), which
     would make the min-detectable-deviation too optimistic.  Instead we
     use the EXACT mean and variance of the Pearson statistic under the
     multinomial alternative (closed form via factorial moments) and
     approximate X^2 ~ g * chi2_h with g, h matched to those moments
     (Patnaik).  This matches simulation to < 0.01 in power.

     Power is a function of the full q, not only delta; to quote a
     1-D power curve in delta we use the representative "shift" family
     q0 = 1/K + eps, qi = (1-q0)/(K-1), for which delta = K*eps/sqrt(K-1).

     min detectable deviation at target power 1-beta:
         delta_min solves power(delta) = 1-beta  (bisection over the
         shift family; closed form only for df=1 under the first-order
         theory: lambda_min = (z_{1-alpha/2} + z_{1-beta})^2)

   Decision rule:
     X^2 > crit                          -> reject H0 (always valid, FPR = alpha)
     X^2 <= crit, power(delta*) >= 1-beta -> accept H0 at resolution delta*
     X^2 <= crit, power(delta*) <  1-beta -> INCONCLUSIVE: report power and
                                             n_needed, not a verdict

   Self-validation below:
     - empirical FPR on a known-uniform RNG ~= alpha
     - analytic power curve ~= empirical power on known-biased generators
     - empirical power at delta_min ~= target power
*/
#include <stdio.h>
#include <math.h>
#include <string.h>

#define K_MAX 64

/* ---------- RNG: xorshift32, known-uniform, deterministic ---------- */
static unsigned rng_state = 2463534242u;
static void rng_seed(unsigned s) { rng_state = s ? s : 2463534242u; }
static unsigned rng_u32(void) {
    unsigned x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_state = x; return x;
}
static double rng_uniform(void) { return (rng_u32() >> 8) * (1.0 / 16777216.0); }

/* ---------- special functions ---------- */
static double log_gamma(double xx) {
    static const double cof[6] = { 76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155,
        0.1208650973866179e-2, -0.5395239384953e-5 };
    double x = xx, y = xx, tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    double ser = 1.000000000190015;
    for (int j = 0; j < 6; j++) { y += 1.0; ser += cof[j] / y; }
    return -tmp + log(2.5066282746310005 * ser / x);
}

/* regularized lower incomplete gamma P(a, x) */
static double gamma_p(double a, double x) {
    if (x <= 0.0) return 0.0;
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int n = 0; n < 1000; n++) {
            ap += 1.0; del *= x / ap; sum += del;
            if (fabs(del) < fabs(sum) * 1e-14) break;
        }
        return sum * exp(-x + a * log(x) - log_gamma(a));
    }
    double b = x + 1.0 - a, c = 1e30, d = 1.0 / b, h = d;
    for (int i = 1; i < 1000; i++) {
        double an = -i * (i - a);
        b += 2.0;
        d = an * d + b; if (fabs(d) < 1e-30) d = 1e-30;
        c = b + an / c; if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double del = d * c; h *= del;
        if (fabs(del - 1.0) < 1e-14) break;
    }
    return 1.0 - exp(-x + a * log(x) - log_gamma(a)) * h;
}

static double chi2_cdf(double x, double df) { return gamma_p(df / 2.0, x / 2.0); }

static double chi2_quantile(double p, double df) {
    double lo = 0.0, hi = df > 1.0 ? df : 1.0;
    while (chi2_cdf(hi, df) < p) hi *= 2.0;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        if (chi2_cdf(mid, df) < p) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* ---------- exact moments of the Pearson statistic under q ----------
   X^2 = sum O_i^2/(N p_i) - N;  multinomial factorial moments give
   E[O_i^2], E[O_i^4], E[O_i^2 O_j^2] in closed form. */
static void exact_moments(const double *p, const double *q, int K, long N,
                          double *mean, double *var) {
    double n = (double)N;
    double n2 = n * (n - 1), n3 = n2 * (n - 2), n4 = n3 * (n - 3);
    double eo2[K_MAX], eo4[K_MAX];
    for (int i = 0; i < K; i++) {
        double qi = q[i];
        eo2[i] = n2 * qi * qi + n * qi;
        eo4[i] = n4 * qi * qi * qi * qi + 6.0 * n3 * qi * qi * qi
               + 7.0 * n2 * qi * qi + n * qi;
    }
    double m = -n, v = 0.0;
    for (int i = 0; i < K; i++) {
        double np = n * p[i];
        m += eo2[i] / np;
        v += (eo4[i] - eo2[i] * eo2[i]) / (np * np);
    }
    for (int i = 0; i < K; i++)
        for (int j = i + 1; j < K; j++) {
            double qi = q[i], qj = q[j];
            double eo2o2 = n4 * qi * qi * qj * qj
                         + n3 * (qi * qi * qj + qi * qj * qj)
                         + n2 * qi * qj;
            double cov = eo2o2 - eo2[i] * eo2[j];
            v += 2.0 * cov / (n * p[i]) / (n * p[j]);
        }
    *mean = m;
    *var = v > 1e-12 ? v : 1e-12;
}

/* Patnaik two-moment approximation: X^2 ~ g * chi2_h */
static double power_vs(const double *p, const double *q, int K, long N,
                       double alpha) {
    double crit = chi2_quantile(1.0 - alpha, K - 1);
    double m, v;
    exact_moments(p, q, K, N, &m, &v);
    double g = v / (2.0 * m), h = 2.0 * m * m / v;
    return 1.0 - chi2_cdf(crit / g, h);
}

/* representative alternative family at divergence delta (uniform shift) */
static void shift_family(double *q, int K, double delta) {
    double eps = delta * sqrt(K - 1.0) / K;
    q[0] = 1.0 / K + eps;
    for (int i = 1; i < K; i++) q[i] = (1.0 - q[0]) / (K - 1);
}

static double gof_power(long N, int K, double alpha, double delta) {
    double p[K_MAX], q[K_MAX];
    for (int i = 0; i < K; i++) p[i] = 1.0 / K;
    shift_family(q, K, delta);
    return power_vs(p, q, K, N, alpha);
}

static double delta_min(long N, int K, double alpha, double target_power) {
    double lo = 0.0, hi = 0.05;
    while (gof_power(N, K, alpha, hi) < target_power) hi *= 2.0;
    for (int i = 0; i < 100; i++) {
        double mid = 0.5 * (lo + hi);
        if (gof_power(N, K, alpha, mid) < target_power) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

static long n_needed(int K, double alpha, double delta, double target_power) {
    long lo = 1, hi = 16;
    while (gof_power(hi, K, alpha, delta) < target_power && hi < (1L << 40)) hi *= 2;
    while (hi - lo > 1) {
        long mid = lo + (hi - lo) / 2;
        if (gof_power(mid, K, alpha, delta) < target_power) lo = mid; else hi = mid;
    }
    return hi;
}

/* ---------- the test ---------- */
typedef struct {
    double stat;        /* chi-square statistic */
    double p_value;
    double crit;
    double power;       /* power at delta_interest given N */
    double delta_min;   /* min detectable delta at target power */
    long   n_needed;    /* N needed for target power at delta_interest */
    int    conclusive;  /* 0 = underpowered: report power, not a verdict */
    int    reject;      /* verdict, meaningful only when conclusive */
} GofResult;

static GofResult gof_test(const long *obs, const double *p, int K, long N,
                          double alpha, double delta_interest, double target_power) {
    GofResult r; memset(&r, 0, sizeof(r));
    double stat = 0.0;
    for (int i = 0; i < K; i++) {
        double e = (double)N * p[i], d = obs[i] - e;
        stat += d * d / e;
    }
    r.stat = stat;
    r.crit = chi2_quantile(1.0 - alpha, K - 1);
    r.p_value = 1.0 - chi2_cdf(stat, K - 1);
    r.power = gof_power(N, K, alpha, delta_interest);
    r.delta_min = delta_min(N, K, alpha, target_power);
    r.n_needed = n_needed(K, alpha, delta_interest, target_power);
    if (stat > r.crit)            { r.conclusive = 1; r.reject = 1; }
    else if (r.power >= target_power) { r.conclusive = 1; r.reject = 0; }
    else                          { r.conclusive = 0; r.reject = 0; }
    return r;
}

/* ---------- known-answer generators for validation ---------- */

/* categorical sampler: bin 0 shifted by eps away from uniform 1/K */
static int sample_biased(double eps, int K) {
    double q0 = 1.0 / K + eps;
    double u = rng_uniform();
    if (u < q0) return 0;
    int b = 1 + (int)((u - q0) / (1.0 - q0) * (K - 1));
    return b > K - 1 ? K - 1 : b;
}

/* chi-square divergence delta of the biased sampler above (exact) */
static double delta_of_eps(double eps, int K) {
    double p = 1.0 / K, q0 = p + eps, qi = (1.0 - q0) / (K - 1);
    return sqrt((q0 - p) * (q0 - p) / p + (K - 1) * (qi - p) * (qi - p) / p);
}

/* empirical rejection rate over M trials of N samples each */
static double empirical_power(long N, int K, double alpha, double eps, int M) {
    double crit = chi2_quantile(1.0 - alpha, K - 1);
    double e = (double)N / K;
    int rej = 0;
    for (int t = 0; t < M; t++) {
        long obs[K_MAX] = { 0 };
        for (long i = 0; i < N; i++) obs[sample_biased(eps, K)]++;
        double stat = 0.0;
        for (int i = 0; i < K; i++) { double d = obs[i] - e; stat += d * d / e; }
        if (stat > crit) rej++;
    }
    return (double)rej / M;
}

/* ---------- application: cosine-weighted hemisphere sampler ---------- */
/* same mapping as material.c; w = cos(theta) has pdf 2w, so t = w^2 is
   uniform on [0,1) under a correct sampler -> uniform bins in t */
static double cosine_t(double u1, int biased) {
    double z = biased ? sqrt(1.0 - u1 * u1) : sqrt(1.0 - u1);
    return z * z;
}

static GofResult test_cosine_sampler(long N, int K, int biased, unsigned seed) {
    rng_seed(seed);
    long obs[K_MAX] = { 0 };
    double p[K_MAX];
    for (int i = 0; i < K; i++) p[i] = 1.0 / K;
    for (long i = 0; i < N; i++) {
        double t = cosine_t(rng_uniform(), biased);
        int b = (int)(t * K);
        obs[b > K - 1 ? K - 1 : b]++;
    }
    return gof_test(obs, p, K, N, 0.05, 0.10, 0.80);
}

/* ---------- self-checks ---------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const int K = 16;
    const double alpha = 0.05;
    const long N = 1000;
    const int M = 4000;

    /* special functions against known values */
    check("cdf_chi2_1", fabs(chi2_cdf(3.8414588, 1) - 0.95) < 1e-4);
    check("cdf_chi2_2", fabs(chi2_cdf(5.9914645, 2) - 0.95) < 1e-4);

    /* exact moments under H0: mean must be exactly K-1 */
    {
        double p[K_MAX], m, v;
        for (int i = 0; i < K; i++) p[i] = 1.0 / K;
        exact_moments(p, p, K, N, &m, &v);
        printf("diag H0 mean=%.6f (K-1=%d) var=%.4f\n", m, K - 1, v);
        check("exact_mean_h0", fabs(m - (K - 1)) < 1e-9);
    }

    /* false positive rate on a known-uniform RNG */
    rng_seed(12345);
    double fpr = empirical_power(N, K, alpha, 0.0, M);
    printf("diag fpr=%.4f alpha=%.2f M=%d\n", fpr, alpha, M);
    check("fpr_controlled", fabs(fpr - alpha) < 0.011); /* 3-sigma window */

    /* analytic power curve vs empirical, at three deviation sizes */
    double eps_pts[3] = { 0.02, 0.03, 0.04 };
    int curve_ok = 1;
    for (int i = 0; i < 3; i++) {
        double delta = delta_of_eps(eps_pts[i], K);
        double analytic = gof_power(N, K, alpha, delta);
        rng_seed(777 + i);
        double emp = empirical_power(N, K, alpha, eps_pts[i], M);
        printf("diag eps=%.3f delta=%.4f analytic=%.4f empirical=%.4f\n",
               eps_pts[i], delta, analytic, emp);
        if (fabs(analytic - emp) > 0.015) curve_ok = 0;
    }
    check("power_curve_matches_empirical", curve_ok);

    /* min detectable deviation: empirical power at delta_min ~= target 0.80 */
    double dm = delta_min(N, K, alpha, 0.80);
    double eps_mdd = dm * sqrt(K - 1.0) / K; /* delta = K*eps/sqrt(K-1) here */
    rng_seed(4242);
    double emp_mdd = empirical_power(N, K, alpha, eps_mdd, M);
    printf("diag delta_min=%.4f eps_equiv=%.4f empirical_power=%.4f target=0.80\n",
           dm, eps_mdd, emp_mdd);
    check("mdd_matches_empirical", fabs(emp_mdd - 0.80) < 0.025);

    /* underpowered: small N must report power, not a verdict */
    {
        rng_seed(99);
        long obs[K_MAX] = { 0 };
        double p[K_MAX];
        for (int i = 0; i < K; i++) p[i] = 1.0 / K;
        long n_small = 100;
        for (long i = 0; i < n_small; i++) obs[(int)(rng_uniform() * K) % K]++;
        GofResult r = gof_test(obs, p, K, n_small, alpha, 0.15, 0.80);
        printf("diag smallN stat=%.2f crit=%.2f power=%.3f n_needed=%ld conclusive=%d\n",
               r.stat, r.crit, r.power, r.n_needed, r.conclusive);
        check("underpowered_reports_power",
              !r.conclusive && r.power < 0.80 && r.n_needed > n_small);
    }

    /* powered and clean: known-uniform source accepted at resolution 0.10 */
    {
        rng_seed(2024);
        long obs[K_MAX] = { 0 };
        double p[K_MAX];
        for (int i = 0; i < K; i++) p[i] = 1.0 / K;
        long n_big = 2500; /* power(0.10) = 0.87 here; 2000 would be borderline */
        for (long i = 0; i < n_big; i++) obs[(int)(rng_uniform() * K) % K]++;
        GofResult r = gof_test(obs, p, K, n_big, alpha, 0.10, 0.80);
        printf("diag bigN stat=%.2f power=%.3f conclusive=%d reject=%d\n",
               r.stat, r.power, r.conclusive, r.reject);
        check("powered_accept_uniform", r.conclusive && !r.reject);
    }

    /* known-biased source rejected */
    {
        rng_seed(31337);
        long obs[K_MAX] = { 0 };
        double p[K_MAX];
        for (int i = 0; i < K; i++) p[i] = 1.0 / K;
        long n_big = 2000;
        for (long i = 0; i < n_big; i++) obs[sample_biased(0.05, K)]++;
        GofResult r = gof_test(obs, p, K, n_big, alpha, 0.10, 0.80);
        printf("diag biased stat=%.2f crit=%.2f reject=%d\n", r.stat, r.crit, r.reject);
        check("known_biased_rejected", r.conclusive && r.reject);
    }

    /* application: cosine-weighted sampler (material.c mapping) */
    {
        GofResult ok = test_cosine_sampler(4000, K, 0, 555);
        printf("diag cosine stat=%.2f crit=%.2f power=%.3f conclusive=%d reject=%d\n",
               ok.stat, ok.crit, ok.power, ok.conclusive, ok.reject);
        check("cosine_sampler_unbiased", ok.conclusive && !ok.reject);

        GofResult bad = test_cosine_sampler(4000, K, 1, 555);
        printf("diag cosine_biased stat=%.2f reject=%d\n", bad.stat, bad.reject);
        check("cosine_biased_detected", bad.conclusive && bad.reject);
    }

    /* power curve table: csv,N,delta,power for K=16, alpha=0.05 */
    long Ns[3] = { 200, 1000, 5000 };
    for (int j = 0; j < 3; j++)
        for (int d = 1; d <= 30; d++)
            printf("csv %ld %.2f %.4f\n", Ns[j], d / 100.0,
                   gof_power(Ns[j], K, alpha, d / 100.0));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
