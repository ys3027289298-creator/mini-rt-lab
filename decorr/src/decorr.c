
/* decorr: per-pixel randomization of a shared low-discrepancy sequence.
 *
 * Each pixel reuses the SAME base 2D Halton point set, but gets independent
 * randomization parameters derived deterministically from (seed, x, y):
 *   - Cranley-Patterson rotation: u'_i = (u_i + xi) mod 1, xi ~ U[0,1)^2
 *   - index permutation:         two independent uniform S_n permutations
 *                                pi_u, pi_v (Fisher-Yates), one per dimension,
 *                                so sample i becomes (u_piu(i), v_piv(i))
 *                                (a shared permutation only reorders points
 *                                 and leaves every average unchanged)
 *
 * Theory (see README.md for the derivation):
 *   The estimator of a pixel is a function of its own randomization parameters.
 *   When two pixels draw independent parameters, their estimators are functions
 *   of independent random variables, hence independent, hence their population
 *   Pearson correlation is exactly 0: rho = 0 (the theoretical upper bound).
 *
 *   On a finite image we measure the sample Pearson correlation r over M
 *   DISJOINT pixel pairs (disjoint => independent pairs under H0). Under H0
 *   E[r] = 0 and Var(r) ~ 1/(M-1), giving distribution-free bounds:
 *     Chebyshev:  P(|r| > 1/sqrt(alpha*(M-1))) <= alpha
 *     3-sigma:    P(|r| > 3/sqrt(M-1))        <= ~0.0027 (normal approximation)
 *
 * With a fixed seed every output is bit-for-bit reproducible.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define N_MAX 256
#define W_MAX 96
#define H_MAX 96

typedef enum { MODE_NONE = 0, MODE_CP = 1, MODE_PERM = 2, MODE_BOTH = 3 } Mode;

typedef struct {
    int width, height, n;
    uint64_t seed;
    int mode; /* OR of MODE_CP / MODE_PERM; MODE_NONE = identical sequence */
} Config;

typedef struct {
    int npix;
    double mean_err;       /* mean signed error vs analytic truth (bias check) */
    double rmse;           /* per-pixel RMSE vs analytic truth */
    double r_h, r_v;       /* sample Pearson corr on disjoint H/V pixel pairs */
    double lowfreq;        /* normalized low-frequency DFT energy of error field */
    double lowfreq_white;  /* expected value for pure white noise: k/N */
} Metrics;

/* ---------- deterministic hashing: all per-pixel randomness comes from here ---------- */

static uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

/* keyed hash, four independent 64-bit channels from (seed, x, y) */
static uint64_t hash64(uint64_t seed, int x, int y, int channel) {
    uint64_t h = seed;
    h = mix64(h ^ (0xd1b54a32d192ed03ULL + (uint64_t)(unsigned)x * 0x9e3779b97f4a7c15ULL));
    h = mix64(h ^ (0x632be597f7f3c2f1ULL + (uint64_t)(unsigned)y * 0xc2b2ae3d27d4eb4fULL));
    h = mix64(h ^ (0x8eb46006c8a11a63ULL + (uint64_t)channel * 0x165667b19e3779f9ULL));
    return h;
}

static double hash_unit(uint64_t seed, int x, int y, int channel) {
    return (double)(hash64(seed, x, y, channel) >> 11) * 0x1.0p-53; /* [0,1) */
}

/* ---------- base sequence: 2D Halton (bases 2 and 3), shared by every pixel ---------- */

static void halton2(unsigned i, double u[2]) {
    double f = 1.0, r = 0.0;
    unsigned k = i;
    while (k > 0) { f /= 2.0; r += f * (k % 2); k /= 2; }
    u[0] = r;

    f = 1.0; r = 0.0; k = i;
    while (k > 0) { f /= 3.0; r += f * (k % 3); k /= 3; }
    u[1] = r;
}

/* ---------- test integrand: non-separable in paired coordinates, analytic truth ----------
 * f_p(u,v) = 1 + a*(u-px)^2 + b*(v-py)^2 + c*(u-px)*(v-py)
 * Integral over [0,1)^2: I = 1 + a*(1/3 - px + px^2) + b*(1/3 - py + py^2)
 *                                       + c*(1/2 - px)*(1/2 - py)
 * The cross term is pairing-sensitive, so a permuted sample index changes the
 * estimate; a purely 1-D average would be permutation-invariant.
 */
static void pixel_params(int x, int y, const Config* cfg, double p[2], double coef[3]) {
    p[0] = (x + 0.5) / cfg->width;
    p[1] = (y + 0.5) / cfg->height;
    coef[0] = 0.6 + 0.3 * p[0];
    coef[1] = 0.5 + 0.3 * p[1];
    coef[2] = 0.7;
}

static double integrand(double u, double v, const double p[2], const double coef[3]) {
    double du = u - p[0], dv = v - p[1];
    return 1.0 + coef[0] * du * du + coef[1] * dv * dv + coef[2] * du * dv;
}

static double truth(const double p[2], const double coef[3]) {
    return 1.0
         + coef[0] * (1.0 / 3.0 - p[0] + p[0] * p[0])
         + coef[1] * (1.0 / 3.0 - p[1] + p[1] * p[1])
         + coef[2] * (0.5 - p[0]) * (0.5 - p[1]);
}

/* ---------- parameter validation: out-of-range values are reported, never clamped ---------- */

static int validate(const Config* cfg) {
    int bad = 0;
    if (cfg->width <= 0 || cfg->width > W_MAX) {
        fprintf(stderr, "param_error: width=%d out of range [1,%d]\n", cfg->width, W_MAX);
        bad = 1;
    }
    if (cfg->height <= 0 || cfg->height > H_MAX) {
        fprintf(stderr, "param_error: height=%d out of range [1,%d]\n", cfg->height, H_MAX);
        bad = 1;
    }
    if (cfg->n <= 1 || cfg->n > N_MAX) {
        fprintf(stderr, "param_error: n=%d out of range [2,%d]\n", cfg->n, N_MAX);
        bad = 1;
    }
    if (cfg->mode < 0 || cfg->mode > MODE_BOTH) {
        fprintf(stderr, "param_error: mode=%d out of range [0,%d]\n", cfg->mode, MODE_BOTH);
        bad = 1;
    }
    return bad;
}

/* ---------- one rendered image + pattern / correlation metrics ---------- */

static double g_err[H_MAX][W_MAX];
static int g_perm_u[N_MAX];
static int g_perm_v[N_MAX];

/* Pearson correlation on DISJOINT pairs: (0,1),(2,3),... => independent pairs */
static double pearson(int hh, int ww, int horizontal) {
    int m2 = horizontal ? (ww / 2) : (hh / 2);
    if (m2 < 2) return 0.0;
    int pairs = horizontal ? hh * m2 : ww * m2;
    double ma = 0.0, mb = 0.0;
    if (horizontal)
        for (int y = 0; y < hh; y++)
            for (int x = 0; x + 1 < ww; x += 2) { ma += g_err[y][x]; mb += g_err[y][x + 1]; }
    else
        for (int y = 0; y + 1 < hh; y += 2)
            for (int x = 0; x < ww; x++) { ma += g_err[y][x]; mb += g_err[y + 1][x]; }
    ma /= pairs; mb /= pairs;
    double cov = 0.0, va = 0.0, vb = 0.0;
    if (horizontal)
        for (int y = 0; y < hh; y++)
            for (int x = 0; x + 1 < ww; x += 2) {
                double a = g_err[y][x] - ma, b = g_err[y][x + 1] - mb;
                cov += a * b; va += a * a; vb += b * b;
            }
    else
        for (int y = 0; y + 1 < hh; y += 2)
            for (int x = 0; x < ww; x++) {
                double a = g_err[y][x] - ma, b = g_err[y + 1][x] - mb;
                cov += a * b; va += a * a; vb += b * b;
            }
    double den = sqrt(va * vb);
    return den > 0.0 ? cov / den : 0.0;
}

static Metrics render(const Config* cfg) {
    Metrics m;
    memset(&m, 0, sizeof(m));

    for (int y = 0; y < cfg->height; y++) {
        for (int x = 0; x < cfg->width; x++) {
            double p[2], coef[3];
            pixel_params(x, y, cfg, p, coef);

            if (cfg->mode & MODE_PERM) {
                /* Two INDEPENDENT uniform S_n permutations, one per dimension:
                 * sample i becomes (u_piu(i), v_piv(i)). A shared index permutation
                 * would merely reorder the same 2D points and leave every average
                 * unchanged; independent per-dimension permutations break the pairing
                 * and therefore change the point set, as scrambling does. */
                for (int i = 0; i < cfg->n; i++) { g_perm_u[i] = i; g_perm_v[i] = i; }
                for (int i = cfg->n - 1; i > 0; i--) {
                    uint64_t hu = hash64(cfg->seed, x, y, 3 + (unsigned)i);
                    uint64_t hv = hash64(cfg->seed, x, y, N_MAX + 3 + (unsigned)i);
                    int ju = (int)(hu % (uint64_t)(i + 1));
                    int jv = (int)(hv % (uint64_t)(i + 1));
                    int tu = g_perm_u[i]; g_perm_u[i] = g_perm_u[ju]; g_perm_u[ju] = tu;
                    int tv = g_perm_v[i]; g_perm_v[i] = g_perm_v[jv]; g_perm_v[jv] = tv;
                }
            }

            double xi[2] = {0.0, 0.0};
            if (cfg->mode & MODE_CP) {
                xi[0] = hash_unit(cfg->seed, x, y, 0);
                xi[1] = hash_unit(cfg->seed, x, y, 1);
            }

            double sum = 0.0;
            for (int i = 0; i < cfg->n; i++) {
                int idx_u = (cfg->mode & MODE_PERM) ? g_perm_u[i] : i;
                int idx_v = (cfg->mode & MODE_PERM) ? g_perm_v[i] : i;
                double u_b[2];
                halton2((unsigned)idx_u, u_b);      /* base-2 component */
                double u_vv;
                double tmp[2];
                halton2((unsigned)idx_v, tmp);      /* base-3 component */
                u_vv = tmp[1];
                double u = u_b[0], v = u_vv;
                if (cfg->mode & MODE_CP) {
                    u += xi[0]; if (u >= 1.0) u -= 1.0;
                    v += xi[1]; if (v >= 1.0) v -= 1.0;
                }
                sum += integrand(u, v, p, coef);
            }
            double est = sum / cfg->n;
            g_err[y][x] = est - truth(p, coef);
        }
    }

    int npix = cfg->width * cfg->height;
    double sum = 0.0, sum2 = 0.0;
    for (int y = 0; y < cfg->height; y++)
        for (int x = 0; x < cfg->width; x++) {
            sum += g_err[y][x]; sum2 += g_err[y][x] * g_err[y][x];
        }
    m.npix = npix;
    m.mean_err = sum / npix;
    m.rmse = sqrt(sum2 / npix);

    m.r_h = pearson(cfg->height, cfg->width, 1);
    m.r_v = pearson(cfg->height, cfg->width, 0);

    /* pattern metric: fraction of DFT energy in the 3x3 lowest-frequency band.
     * Parseval: sum |E(k)|^2 / N^2 = mean e^2, so we only need 9 coefficients.
     * For white noise every band carries k/N of the energy (9/N here). */
    double low_energy = 0.0;
    for (int ky = -1; ky <= 1; ky++)
        for (int kx = -1; kx <= 1; kx++) {
            double re = 0.0, im = 0.0;
            for (int y = 0; y < cfg->height; y++)
                for (int x = 0; x < cfg->width; x++) {
                    double ang = 2.0 * PI * ((double)kx * x / cfg->width
                                           + (double)ky * y / cfg->height);
                    re += g_err[y][x] * cos(ang);
                    im -= g_err[y][x] * sin(ang);
                }
            low_energy += re * re + im * im;
        }
    double total = sum2 * npix; /* sum |DFT|^2 = N * sum e^2 */
    m.lowfreq = total > 0.0 ? low_energy / total : 0.0;
    m.lowfreq_white = 9.0 / npix;
    return m;
}

/* finite-sample theoretical bounds for |r| given M independent disjoint pairs */
static double bound_chebyshev(int pairs, double alpha) {
    return 1.0 / sqrt(alpha * (pairs - 1));
}
static double bound_3sigma(int pairs) {
    return 3.0 / sqrt((double)pairs - 1.0);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const int W = 64, H = 32, N = 24;
    int pairs_h = H * (W / 2);          /* 1024 independent horizontal pairs */
    int pairs_v = W * (H / 2);          /* 1024 independent vertical pairs   */
    double t_h_cheb = bound_chebyshev(pairs_h, 0.01);
    double t_h_3sig = bound_3sigma(pairs_h);
    double t_v_cheb = bound_chebyshev(pairs_v, 0.01);

    /* --- 1. parameter validation: illegal values must be reported, not clamped --- */
    Config bad_cfg[] = {
        {0, 32, 24, 1, MODE_CP},        /* width = 0   */
        {64, -1, 24, 1, MODE_CP},       /* height < 0  */
        {64, 32, 1, 1, MODE_CP},        /* n = 1       */
        {64, 32, 300, 1, MODE_CP},      /* n > N_MAX   */
        {64, 32, 24, 1, 7},             /* illegal mode */
    };
    int nbad = (int)(sizeof(bad_cfg) / sizeof(bad_cfg[0]));
    int rejected = 0;
    for (int i = 0; i < nbad; i++)
        rejected += validate(&bad_cfg[i]);
    check("out_of_range_reported", rejected == nbad);
    Config good = {W, H, N, 1, MODE_CP};
    check("valid_config_accepted", validate(&good) == 0);

    /* --- 2. headline comparison across modes (fixed seed => deterministic) --- */
    const char* names[4] = {"none", "cp", "perm", "cp+perm"};
    Metrics mm[4];
    printf("\n%-8s %10s %10s %10s %10s %12s %12s\n",
           "mode", "mean_err", "rmse", "r_h", "r_v", "lowfreq", "white_ref");
    for (int mode = 0; mode <= 3; mode++) {
        Config cfg = {W, H, N, 12345, mode};
        mm[mode] = render(&cfg);
        printf("%-8s %10.2e %10.2e %10.4f %10.4f %12.5f %12.5f\n",
               names[mode], mm[mode].mean_err, mm[mode].rmse,
               mm[mode].r_h, mm[mode].r_v, mm[mode].lowfreq, mm[mode].lowfreq_white);
    }

    /* --- 3. theory bound: rho = 0 for independent randomization parameters,
     *     measured |r| must stay under the finite-sample bounds.
     *     CP rotation fully breaks the shared point structure; permutation alone
     *     only reorders fixed marginals, so its decorrelation is partial (its
     *     estimate still contains the common marginal quadrature bias). --- */
    printf("\nbound_chebyshev_99=%.4f bound_3sigma=%.4f (M=%d disjoint pairs, rho=0)\n",
           t_h_cheb, t_h_3sig, pairs_h);
    check("corr_none_above_bound", fabs(mm[MODE_NONE].r_h) > t_h_cheb); /* fixed pattern: correlated */
    for (int mode = 1; mode <= 3; mode++) {
        if (!(mode & MODE_CP)) continue; /* full rho=0 claim needs the rotation */
        char label[64];
        snprintf(label, sizeof(label), "corr_%s_under_chebyshev", names[mode]);
        check(label, fabs(mm[mode].r_h) <= t_h_cheb && fabs(mm[mode].r_v) <= t_v_cheb);
    }
    /* permutation alone must still decorrelate strictly, just not all the way */
    check("corr_perm_partial",
          fabs(mm[MODE_PERM].r_h) < fabs(mm[MODE_NONE].r_h) - 0.05);

    /* --- 4. pattern metric must drop relative to the shared-sequence baseline --- */
    for (int mode = 1; mode <= 3; mode++) {
        if (mode == MODE_PERM) continue; /* marginals preserve most low-freq energy */
        char label[64];
        snprintf(label, sizeof(label), "pattern_%s_drops", names[mode]);
        check(label, mm[mode].lowfreq < 0.5 * mm[MODE_NONE].lowfreq);
    }
    check("pattern_perm_drops",
          mm[MODE_PERM].lowfreq < mm[MODE_NONE].lowfreq);
    /* the full combination reaches the white-noise spectral floor */
    check("pattern_near_white_floor", mm[MODE_BOTH].lowfreq < 3.0 * mm[MODE_BOTH].lowfreq_white);

    /* --- 5. CP rotation is unbiased on a single image (random shift expectation
     *     preserves the integral); perm-only is not, because its marginal point
     *     sets are fixed and carry quadrature bias regardless of pairing. --- */
    for (int mode = 1; mode <= 3; mode++) {
        if (!(mode & MODE_CP)) continue;
        char label[64];
        snprintf(label, sizeof(label), "unbiased_%s", names[mode]);
        double tol = 3.0 * mm[mode].rmse / sqrt((double)mm[mode].npix);
        check(label, fabs(mm[mode].mean_err) <= tol);
    }

    /* --- 6. reproducibility: fixed params must produce identical results --- */
    Config cfg_a = {W, H, N, 999, MODE_BOTH};
    Config cfg_b = cfg_a;
    Metrics a = render(&cfg_a), b = render(&cfg_b);
    check("reproducible_same_seed",
          a.mean_err == b.mean_err && a.rmse == b.rmse &&
          a.r_h == b.r_h && a.r_v == b.r_v && a.lowfreq == b.lowfreq);
    Config cfg_c = cfg_a; cfg_c.seed = 1000;
    Metrics c = render(&cfg_c);
    check("seed_changes_output",
          c.r_h != a.r_h || c.lowfreq != a.lowfreq);

    /* --- 7. multi-seed trials: violation rate must match the bound's alpha level.
     *     Only modes containing CP rotation satisfy the full rho=0 premise;
     *     perm-only keeps a shared marginal-bias component by construction. --- */
    const int TRIALS = 32;
    double alpha = 0.05;
    double t_cheb05 = bound_chebyshev(pairs_h, alpha);
    const int test_modes[2] = {MODE_CP, MODE_BOTH};
    int viol_cheb[2] = {0, 0}, viol_3sig[2] = {0, 0};
    double abs_r_avg[2] = {0, 0};
    for (int t = 0; t < TRIALS; t++) {
        for (int k = 0; k < 2; k++) {
            Config cfg = {W, H, N, 7000 + (uint64_t)t * 131, test_modes[k]};
            Metrics r = render(&cfg);
            double ar = fabs(r.r_h);
            abs_r_avg[k] += ar;
            if (ar > t_cheb05) viol_cheb[k]++;
            if (ar > t_h_3sig) viol_3sig[k]++;
        }
    }
    printf("\ntrials=%d bound_chebyshev_95=%.4f bound_3sigma=%.4f\n",
           TRIALS, t_cheb05, t_h_3sig);
    int rate_ok = 1;
    for (int k = 0; k < 2; k++) {
        abs_r_avg[k] /= TRIALS;
        printf("%-8s mean_|r|=%.4f chebyshev_viol=%d/%d threesigma_viol=%d/%d\n",
               names[test_modes[k]], abs_r_avg[k], viol_cheb[k], TRIALS,
               viol_3sig[k], TRIALS);
        /* Chebyshev promises <= 5%; allow binomial slack (mean |r| ~ 1/sqrt(M*pi/2)) */
        rate_ok = rate_ok && (double)viol_cheb[k] / TRIALS <= 0.25
                           && abs_r_avg[k] < 0.05;
    }
    check("violation_rate_within_bound", rate_ok);

    /* --- 7b. what permutation does and does not remove. A per-dimension permutation
     *      leaves EACH marginal point set unchanged (the same N Halton values per
     *      axis), so additive (marginal) quadrature bias is preserved under every
     *      permutation and cannot be averaged away; only the paired cross-structure
     *      is randomized. Verify directly: over an ensemble of perm images, the
     *      mean error stays near the single-image value (shared marginal bias),
     *      while CP images ensemble-average to ~0 because rotation moves marginals. */
    double perm_ens = 0.0, cp_ens = 0.0;
    const int ENS = 64;
    for (int t = 0; t < ENS; t++) {
        Config cp_cfg = {W, H, N, 20000 + (uint64_t)t * 977, MODE_CP};
        Config pm_cfg = {W, H, N, 20000 + (uint64_t)t * 977, MODE_PERM};
        cp_ens += render(&cp_cfg).mean_err;
        perm_ens += render(&pm_cfg).mean_err;
    }
    cp_ens /= ENS; perm_ens /= ENS;
    printf("ensemble_mean_err: perm=%.3e (shared marginal bias, single=%.3e)  cp=%.3e\n",
           perm_ens, mm[MODE_PERM].mean_err, cp_ens);
    check("perm_keeps_marginal_bias",
          fabs(perm_ens - mm[MODE_PERM].mean_err) < 0.25 * fabs(mm[MODE_PERM].mean_err));
    check("cp_ensemble_unbiased", fabs(cp_ens) < 0.25 * fabs(mm[MODE_PERM].mean_err));

    /* --- 8. hash sanity: CP shifts of one pixel across seeds must be uniform on [0,1) --- */
    int bins = 16, cnt[16] = {0};
    int draws = 16000;
    for (int t = 0; t < draws; t++) {
        double u = hash_unit((uint64_t)t, 7, 9, 0);
        int b = (int)(u * bins);
        if (b >= 0 && b < bins) cnt[b]++;
    }
    double chi = 0.0, expected = (double)draws / bins;
    for (int i = 0; i < bins; i++)
        chi += (cnt[i] - expected) * (cnt[i] - expected) / expected;
    /* df=15, 99% critical value 30.58 */
    printf("hash_uniformity_chi2=%.2f (df=15, crit99=30.58)\n", chi);
    check("hash_shift_uniform", chi < 30.58);

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
