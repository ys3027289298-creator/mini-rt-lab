/* photon: bandwidth selection for photon-mapping density estimation (KDE).
 *
 * Radiance estimate L(x) ~ (1/(n h^d)) sum K((x-x_p)/h) Phi_p is a kernel
 * density estimate. With a second-order kernel:
 *   bias(x)     ~ (h^2 / 2) mu2(K) * laplacian f(x)     -> O(h^2)
 *   variance(x) ~ R(K) f(x) / (n h^d)                   -> O(1/(n h^d))
 * Balancing bias^2 against variance gives the AMISE-optimal rate
 * h* ~ n^(-1/(d+4)); the Silverman rule supplies the constant:
 *   h = sigma * (4/(d+2))^(1/(d+4)) * n^(-1/(d+4))
 * The self-test checks these relations, measures the error against an
 * analytic density, and compares the rule-based bandwidth with an
 * exhaustive search over h.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

static const double PI = 3.14159265358979323846;

/* standard normal kernel moments: mu2 = integral u^2 K(u) du,
   RK = integral K(u)^2 du = (4 pi)^(-d/2) for the product kernel */
static double kernel_mu2(void) { return 1.0; }
static double kernel_RK(int dim) { return pow(4.0 * PI, -(double)dim / 2.0); }

/* Silverman bandwidth from photon count; returns < 0 on invalid input. */
double photon_bandwidth(int n, int dim, double sigma) {
    if (n <= 0 || dim < 1 || !(sigma > 0.0)) return -1.0;
    double e = 1.0 / ((double)dim + 4.0);
    return sigma * pow(4.0 / ((double)dim + 2.0), e) * pow((double)n, -e);
}

/* leading-order bias magnitude upper bound: (h^2/2) mu2 |laplacian f|.
   Tends to 0 as h -> 0; it cannot diverge. */
double photon_bias_bound(double h, double laplacian) {
    if (!(h >= 0.0)) return -1.0;
    return 0.5 * h * h * kernel_mu2() * fabs(laplacian);
}

/* leading-order variance: R(K) f / (n h^d); returns < 0 on invalid input. */
double photon_variance(double h, int n, int dim, double f) {
    if (n <= 0 || !(h > 0.0) || dim < 1) return -1.0;
    return kernel_RK(dim) * f / ((double)n * pow(h, (double)dim));
}

/* MSE estimate = bias^2 + variance. */
double photon_mse(double h, int n, int dim, double f, double laplacian) {
    double b = photon_bias_bound(h, laplacian);
    double v = photon_variance(h, n, dim, f);
    if (b < 0.0 || v < 0.0) return -1.0;
    return b * b + v;
}

/* ---- self-test ---- */

static double halton(unsigned i, unsigned base) {
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= base; r += f * (double)(i % base); i /= base; }
    return r;
}

static double normal_pdf(double x) { return exp(-0.5 * x * x) / sqrt(2.0 * PI); }

/* Box-Muller on Halton(2,3): deterministic standard-normal samples. */
static double normal_sample(unsigned i) {
    double u1 = halton(i + 1, 2);
    double u2 = halton(i + 1, 3);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

static double kde(const double* xs, int n, double h, double x) {
    double acc = 0.0;
    for (int i = 0; i < n; i++) {
        double u = (x - xs[i]) / h;
        acc += exp(-0.5 * u * u);
    }
    return acc / ((double)n * h * sqrt(2.0 * PI));
}

/* RMSE of the KDE against the analytic N(0,1) density on a grid. */
static double kde_rmse(const double* xs, int n, double h) {
    double acc = 0.0;
    int m = 0;
    for (int j = 0; j <= 120; j++) {
        double x = -3.0 + 6.0 * (double)j / 120.0;
        double d = kde(xs, n, h, x) - normal_pdf(x);
        acc += d * d;
        m++;
    }
    return sqrt(acc / (double)m);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* zero or negative photon count must be rejected */
    check("reject_zero_photons", photon_bandwidth(0, 2, 1.0) < 0.0);
    check("reject_negative_photons", photon_bandwidth(-16, 2, 1.0) < 0.0);
    check("reject_zero_variance", photon_variance(0.1, 0, 2, 1.0) < 0.0);
    check("reject_zero_h", photon_variance(0.0, 64, 2, 1.0) < 0.0);

    /* bandwidth from photon count: h ~ n^(-1/(d+4)) */
    double h1 = photon_bandwidth(1000, 2, 1.0);
    double h8 = photon_bandwidth(8000, 2, 1.0);
    check("bandwidth_positive", h1 > 0.0 && h8 > 0.0);
    check("bandwidth_shrinks", h8 < h1);
    check("bandwidth_rate", fabs(h8 / h1 - pow(8.0, -1.0 / 6.0)) < 1e-9);

    /* bias ~ h^2, variance ~ 1/(n h^d) */
    double b_h = photon_bias_bound(0.3, 1.0);
    double b_2h = photon_bias_bound(0.6, 1.0);
    check("bias_quadratic", fabs(b_2h / b_h - 4.0) < 1e-12);
    double v_n = photon_variance(0.3, 100, 2, 1.0);
    double v_2n = photon_variance(0.3, 200, 2, 1.0);
    double v_2h = photon_variance(0.6, 100, 2, 1.0);
    check("var_inverse_n", fabs(v_2n / v_n - 0.5) < 1e-12);
    check("var_inverse_hd", fabs(v_2h / v_n - 0.25) < 1e-12);

    /* bias bound must not diverge as h -> 0: strictly decreasing to 0 */
    double prev = 1e30;
    int monotone = 1;
    double h = 0.5;
    for (int k = 0; k < 24; k++) {
        double b = photon_bias_bound(h, 2.0);
        if (!(b < prev) || !isfinite(b)) monotone = 0;
        prev = b;
        h *= 0.5;
    }
    check("bias_vanishes_h_to_zero", monotone && prev < 1e-9);

    /* bias/variance trade-off table (printed for inspection) */
    printf("# h bias_bound std\n");
    for (int k = 0; k < 5; k++) {
        double hh = 0.1 * (double)(1 << k);
        printf("# %.3f %.6f %.6f\n", hh, photon_bias_bound(hh, 0.4),
               sqrt(photon_variance(hh, 1024, 1, 0.4)));
    }

    /* empirical error vs analytic N(0,1) density */
    enum { N = 1024, RUNS = 16, NSUB = 256, HSTEPS = 48 };
    static double xs[N];
    for (int i = 0; i < N; i++) xs[i] = normal_sample((unsigned)i);

    double h_rule = photon_bandwidth(N, 1, 1.0);
    double err_rule = kde_rmse(xs, N, h_rule);

    /* exhaustive search for the best bandwidth on a log grid */
    double best_err = 1e30, best_h = 0.0;
    for (int s = 0; s < HSTEPS; s++) {
        double hh = 0.05 * pow(1.0 / 0.05, (double)s / (double)(HSTEPS - 1));
        double e = kde_rmse(xs, N, hh);
        if (e < best_err) { best_err = e; best_h = hh; }
    }
    printf("# h_rule=%.4f err_rule=%.5f h_best=%.4f err_best=%.5f\n",
           h_rule, err_rule, best_h, best_err);
    check("rule_near_exhaustive", err_rule <= 3.0 * best_err);

    /* measured bias at x=0, averaged over independent runs, must stay
       under the analytic bound (with slack for higher-order terms) */
    double h_sub = photon_bandwidth(NSUB, 1, 1.0);
    double acc = 0.0;
    for (int r = 0; r < RUNS; r++) {
        double sub[NSUB];
        for (int i = 0; i < NSUB; i++)
            sub[i] = normal_sample((unsigned)(r * NSUB + i));
        acc += kde(sub, NSUB, h_sub, 0.0);
    }
    double bias_meas = fabs(acc / RUNS - normal_pdf(0.0));
    /* f''(0) = -phi(0) for the standard normal */
    double bias_ub = photon_bias_bound(h_sub, normal_pdf(0.0));
    printf("# bias_meas=%.5f bias_bound=%.5f\n", bias_meas, bias_ub);
    check("bias_under_bound", bias_meas <= 2.0 * bias_ub);

    /* MSE model has an interior minimum near the rule bandwidth */
    double mse_rule = photon_mse(h_rule, N, 1, normal_pdf(0.0), normal_pdf(0.0));
    double mse_tiny = photon_mse(h_rule * 0.05, N, 1, normal_pdf(0.0), normal_pdf(0.0));
    double mse_huge = photon_mse(h_rule * 8.0, N, 1, normal_pdf(0.0), normal_pdf(0.0));
    check("mse_interior_min", mse_rule < mse_tiny && mse_rule < mse_huge);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
