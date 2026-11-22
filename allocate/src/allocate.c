
/* allocate: per-pixel sample budget allocation from pilot variance.
 *
 * Uniform samples-per-pixel starves noisy pixels and wastes samples on smooth
 * ones. Given per-pixel variance from a small pilot, this module:
 *   - pre-allocates a fixed total budget by Neyman allocation (n_i ~ sigma_i)
 *     with a per-pixel floor and exact integer accounting;
 *   - solves the budget required to reach a target RMSE and drives a
 *     deterministic convergence loop (fixed-seed LCG, documented stop rule);
 *   - predicts adaptive-vs-uniform RMSE and their ratio at equal budget;
 *   - rejects budgets below the minimum usable value and reports that value.
 */
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#define NPIX_MAX 64
#define MAX_ROUNDS 16
#define TWO_PI 6.28318530717958647692

/* return codes */
enum { ALLOC_OK = 0, ALLOC_E_BUDGET = 1, ALLOC_E_TARGET = 2, ALLOC_E_PARAM = 3 };

/* ---------------- deterministic PRNG: fixed-seed LCG, no rand()/time() ---- */
typedef struct { uint32_t s; } Rng;

static void rng_seed(Rng* r, uint32_t seed) { r->s = seed ? seed : 1u; }

static double rng_next(Rng* r) {
    r->s = r->s * 1664525u + 1013904223u;      /* Numerical Recipes LCG */
    return (double)(r->s >> 8) * (1.0 / 16777216.0); /* [0,1), 24 bits */
}

static double rng_normal(Rng* r) {             /* Box-Muller, N(0,1) */
    double u1 = 1.0 - rng_next(r);
    double u2 = rng_next(r);
    return sqrt(-2.0 * log(u1)) * cos(TWO_PI * u2);
}

/* ---------------- synthetic pixel: stands in for tracing one camera ray ---- */
typedef struct { double mu; double sigma; } PixelModel;

static double pixel_sample(const PixelModel* p, Rng* r) {
    return p->mu + p->sigma * rng_normal(r);
}

/* ---------------- Welford online mean / unbiased variance ----------------- */
typedef struct { double mean; double m2; int n; } WStat;

static void stat_push(WStat* s, double x) {
    s->n++;
    double d = x - s->mean;
    s->mean += d / (double)s->n;
    s->m2 += d * (x - s->mean);
}

static double stat_var(const WStat* s) {
    return s->n > 1 ? s->m2 / (double)(s->n - 1) : 0.0;
}

/* ---------------- minimum usable budget ----------------------------------- */
/* Variance needs >= 2 samples per pixel; the caller chooses n_min >= 2. */
static long min_budget(int npix, int n_min) {
    return (long)npix * (long)n_min;
}

/* ---------------- Neyman allocation with per-pixel floor ------------------ */
/* n_i proportional to sigma_i. Pixels whose share would drop below n_min are
 * pinned to n_min and their weight is redistributed (water-filling). The full
 * budget is always spent exactly; leftover units go to the largest fractional
 * remainder, ties to the smallest current count, then the lowest index. */
static int allocate_neyman(int npix, const double sigma[], long budget, int n_min,
                           long out_n[], double out_frac[], long* min_required) {
    long bmin = min_budget(npix, n_min);
    if (min_required) *min_required = bmin;
    if (npix <= 0 || npix > NPIX_MAX || n_min < 2) return ALLOC_E_PARAM;
    if (budget < bmin) return ALLOC_E_BUDGET;           /* zero lands here */

    int active[NPIX_MAX];
    double w[NPIX_MAX], rem[NPIX_MAX];
    for (int i = 0; i < npix; i++) {
        active[i] = 1;
        w[i] = sigma[i] > 0.0 ? sigma[i] : 0.0;
        out_n[i] = 0;
        rem[i] = -1.0;
    }

    long remaining = budget;
    for (;;) {
        double sw = 0.0;
        for (int i = 0; i < npix; i++)
            if (active[i]) sw += w[i];

        if (sw <= 0.0) { /* all active pixels have zero estimated variance */
            for (int i = 0; i < npix; i++)
                if (active[i]) { out_n[i] = n_min; remaining -= n_min; active[i] = 0; }
            break;
        }

        double q[NPIX_MAX];
        int any_below = 0;
        for (int i = 0; i < npix; i++) {
            if (!active[i]) continue;
            q[i] = (double)remaining * w[i] / sw;
            if (q[i] < (double)n_min) any_below = 1;
        }
        if (!any_below) {
            for (int i = 0; i < npix; i++) {
                if (!active[i]) continue;
                double fl = floor(q[i]);
                out_n[i] = (long)fl;
                rem[i] = q[i] - fl;
            }
            break;
        }
        for (int i = 0; i < npix; i++) {
            if (active[i] && q[i] < (double)n_min) {
                out_n[i] = n_min;
                active[i] = 0;
                remaining -= n_min;
            }
        }
    }

    long used = 0;
    for (int i = 0; i < npix; i++) used += out_n[i];
    long leftover = budget - used;
    while (leftover > 0) {
        int best = 0;
        for (int i = 1; i < npix; i++) {
            if (rem[i] > rem[best] + 1e-15) { best = i; continue; }
            if (fabs(rem[i] - rem[best]) <= 1e-15 &&
                (out_n[i] < out_n[best] ||
                 (out_n[i] == out_n[best] && i < best)))
                best = i;
        }
        out_n[best]++;
        rem[best] -= 1.0; /* consumed its remainder: strict largest-remainder order */
        leftover--;
    }

    for (int i = 0; i < npix; i++)
        out_frac[i] = (double)out_n[i] / (double)budget;
    return ALLOC_OK;
}

/* ---------------- predicted error ----------------------------------------- */
/* RMSE of the image estimate when pixel i gets n_i samples with noise sigma_i. */
static double predict_rmse(int npix, const double sigma[], const long n[]) {
    double acc = 0.0;
    for (int i = 0; i < npix; i++)
        acc += sigma[i] * sigma[i] / (double)n[i];
    return sqrt(acc / (double)npix);
}

/* Uniform split n_i = B/N: RMSE_u^2 = sum(sigma_i^2) / B. */
static double predict_uniform_rmse(int npix, const double sigma[], long budget) {
    double acc = 0.0;
    for (int i = 0; i < npix; i++) acc += sigma[i] * sigma[i];
    return sqrt(acc / (double)budget);
}

/* Neyman allocation makes RMSE^2 = S^2 / (N B), so the budget hitting target
 * error eps is B* = S^2 / (N eps^2), clamped up to the per-pixel floor. */
static long budget_for_target(int npix, const double sigma[],
                              double target_rmse, int n_min) {
    if (target_rmse <= 0.0) return -1;
    double s = 0.0;
    for (int i = 0; i < npix; i++)
        if (sigma[i] > 0.0) s += sigma[i];
    double b = s * s / ((double)npix * target_rmse * target_rmse);
    long need = (long)ceil(b);
    long bmin = min_budget(npix, n_min);
    return need < bmin ? bmin : need;
}

/* ---------------- convergence driver -------------------------------------- */
typedef struct {
    int npix;
    int n_min;
    double target_rmse;
    long max_budget;
    uint32_t seed;
} ConvConfig;

typedef struct {
    WStat st[NPIX_MAX];
    long n[NPIX_MAX];
    double frac[NPIX_MAX];
    long budget;
    long min_budget;
    long required_budget; /* populated on ALLOC_E_TARGET */
    double rmse;
    int rounds;
    int converged;
} ConvResult;

/* Stop rule (deterministic, reruns are bit-identical for the same seed):
 *  1. pilot: n_min samples per pixel;
 *  2. each round re-estimates s_i^2 (Welford), runs Neyman at the current
 *     budget, and traces only the delta (sample counts never shrink);
 *  3. predicted RMSE = sqrt(mean_i s_i^2 / n_i);
 *  4. converged when predicted RMSE <= target for two consecutive rounds;
 *  5. budget grows to max(2B, B*) until max_budget; hard stop at MAX_ROUNDS. */
static int converge(const PixelModel* scene, const ConvConfig* cfg, ConvResult* res) {
    memset(res, 0, sizeof(*res));
    if (cfg->npix <= 0 || cfg->npix > NPIX_MAX ||
        cfg->n_min < 2 || cfg->target_rmse <= 0.0)
        return ALLOC_E_PARAM;

    res->min_budget = min_budget(cfg->npix, cfg->n_min);
    if (cfg->max_budget < res->min_budget) {
        res->required_budget = res->min_budget;
        return ALLOC_E_BUDGET;
    }

    Rng r;
    rng_seed(&r, cfg->seed);
    double sigma[NPIX_MAX];
    for (int i = 0; i < cfg->npix; i++)
        for (int k = 0; k < cfg->n_min; k++)
            stat_push(&res->st[i], pixel_sample(&scene[i], &r));

    long budget = res->min_budget;
    int hold = 0;
    for (int round = 1; round <= MAX_ROUNDS; round++) {
        for (int i = 0; i < cfg->npix; i++)
            sigma[i] = sqrt(fmax(stat_var(&res->st[i]), 0.0));

        int rc = allocate_neyman(cfg->npix, sigma, budget, cfg->n_min,
                                 res->n, res->frac, &res->min_budget);
        if (rc != ALLOC_OK) return rc;

        for (int i = 0; i < cfg->npix; i++)
            while ((long)res->st[i].n < res->n[i])
                stat_push(&res->st[i], pixel_sample(&scene[i], &r));

        for (int i = 0; i < cfg->npix; i++)
            sigma[i] = sqrt(fmax(stat_var(&res->st[i]), 0.0));

        res->rmse = predict_rmse(cfg->npix, sigma, res->n);
        res->budget = budget;
        res->rounds = round;

        if (res->rmse <= cfg->target_rmse) {
            if (++hold >= 2) { res->converged = 1; return ALLOC_OK; }
        } else {
            hold = 0;
        }

        long need = budget_for_target(cfg->npix, sigma,
                                      cfg->target_rmse, cfg->n_min);
        res->required_budget = need;
        if (need > cfg->max_budget) return ALLOC_E_TARGET;

        long next = budget * 2;
        if (next < need) next = need;
        if (next > cfg->max_budget) next = cfg->max_budget;
        if (next <= budget) {
            if (hold > 0) continue; /* re-check at cap: identical estimate => hold reaches 2 */
            return ALLOC_E_TARGET;
        }
        budget = next;
    }

    for (int i = 0; i < cfg->npix; i++)
        sigma[i] = sqrt(fmax(stat_var(&res->st[i]), 0.0));
    res->required_budget = budget_for_target(cfg->npix, sigma,
                                             cfg->target_rmse, cfg->n_min);
    return ALLOC_E_TARGET;
}

/* ---------------- empirical end-to-end check (fixed seeds) ---------------- */
static void empirical_mse(const PixelModel* scene, int npix, const long n_adapt[],
                          long n_uniform, int trials, uint32_t seed,
                          double* mse_adapt, double* mse_uniform) {
    Rng ra, ru;
    rng_seed(&ra, seed);
    rng_seed(&ru, seed ^ 0x9e3779b9u);
    double sa = 0.0, su = 0.0;
    for (int t = 0; t < trials; t++) {
        for (int i = 0; i < npix; i++) {
            double ma = 0.0, mu = 0.0;
            for (long k = 0; k < n_adapt[i]; k++) ma += pixel_sample(&scene[i], &ra);
            for (long k = 0; k < n_uniform; k++) mu += pixel_sample(&scene[i], &ru);
            ma /= (double)n_adapt[i];
            mu /= (double)n_uniform;
            double da = ma - scene[i].mu;
            double du = mu - scene[i].mu;
            sa += da * da;
            su += du * du;
        }
    }
    *mse_adapt = sa / (double)(trials * npix);
    *mse_uniform = su / (double)(trials * npix);
}

/* ---------------- self-test ------------------------------------------------ */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    enum { N = 8, NMIN = 4 };
    static const double SIGMA[N] = { 0.5, 0.5, 1.0, 1.0, 2.0, 2.0, 4.0, 4.0 };
    PixelModel scene[N];
    for (int i = 0; i < N; i++) { scene[i].mu = (double)(i + 1); scene[i].sigma = SIGMA[i]; }

    long n[N];
    double frac[N];
    long minreq = -1;

    /* zero / negative / below-minimum budgets are rejected, with the minimum
     * usable value reported */
    check("zero_budget_rejected",
          allocate_neyman(N, SIGMA, 0, NMIN, n, frac, &minreq) == ALLOC_E_BUDGET &&
          minreq == (long)N * NMIN);
    check("negative_budget_rejected",
          allocate_neyman(N, SIGMA, -5, NMIN, n, frac, &minreq) == ALLOC_E_BUDGET &&
          minreq == (long)N * NMIN);
    check("below_min_rejected",
          allocate_neyman(N, SIGMA, (long)N * NMIN - 1, NMIN, n, frac, &minreq)
          == ALLOC_E_BUDGET);
    check("min_budget_ok",
          allocate_neyman(N, SIGMA, (long)N * NMIN, NMIN, n, frac, &minreq) == ALLOC_OK);

    /* exact accounting: whole budget spent, floor respected, fractions sum to 1 */
    long budgets[] = { (long)N * NMIN, (long)N * NMIN + 1, (long)N * NMIN + 7,
                       512, 513, 1000, 4096 };
    int accounting_ok = 1, floor_ok = 1, frac_ok = 1;
    for (unsigned b = 0; b < sizeof(budgets) / sizeof(budgets[0]); b++) {
        long B = budgets[b];
        if (allocate_neyman(N, SIGMA, B, NMIN, n, frac, &minreq) != ALLOC_OK) {
            accounting_ok = 0; continue;
        }
        long sum = 0;
        double fs = 0.0;
        for (int i = 0; i < N; i++) { sum += n[i]; fs += frac[i]; if (n[i] < NMIN) floor_ok = 0; }
        if (sum != B) accounting_ok = 0;
        if (fabs(fs - 1.0) > 1e-12) frac_ok = 0;
    }
    check("alloc_sums_to_budget", accounting_ok);
    check("floor_respected", floor_ok);
    check("fractions_sum_one", frac_ok);

    /* Neyman ratio: sigma 3:1 should buy ~3:1 samples */
    double s2[2] = { 3.0, 1.0 };
    long n2[2];
    double f2[2];
    allocate_neyman(2, s2, 400, 2, n2, f2, &minreq);
    double ratio = (double)n2[0] / (double)n2[1];
    check("neyman_ratio_3to1", ratio > 2.9 && ratio < 3.1);

    /* zero-variance pixel gets the floor and does not break accounting */
    double sz[3] = { 0.0, 1.0, 3.0 };
    long nz[3];
    double fz[3];
    check("zero_sigma_handled",
          allocate_neyman(3, sz, 60, NMIN, nz, fz, &minreq) == ALLOC_OK &&
          nz[0] == NMIN && (nz[0] + nz[1] + nz[2]) == 60);

    /* predicted error: adaptive vs uniform at equal budget */
    long B = 512;
    allocate_neyman(N, SIGMA, B, NMIN, n, frac, &minreq);
    double rmse_a = predict_rmse(N, SIGMA, n);
    double rmse_u = predict_uniform_rmse(N, SIGMA, B);
    double err_ratio = rmse_a / rmse_u;
    check("error_ratio_le_one", err_ratio <= 1.0 + 1e-12);
    check("error_ratio_scene_gain", err_ratio < 0.9);

    /* budget solved for a target actually reaches it */
    double target = 0.5;
    long need = budget_for_target(N, SIGMA, target, NMIN);
    long nt[N];
    double ft[N];
    int tbc = allocate_neyman(N, SIGMA, need, NMIN, nt, ft, &minreq);
    check("target_budget_reaches",
          tbc == ALLOC_OK && predict_rmse(N, SIGMA, nt) <= target);

    /* unreachable target under a tight cap is refused with the required budget */
    ConvConfig bad = { N, NMIN, 1e-3, (long)N * NMIN, 12345u };
    ConvResult rb;
    check("target_unreachable_rejected",
          converge(scene, &bad, &rb) == ALLOC_E_TARGET &&
          rb.required_budget > bad.max_budget);

    /* convergence succeeds at an achievable target */
    ConvConfig cfg = { N, NMIN, 0.8, (long)N * 2048, 12345u };
    ConvResult r1, r2;
    int rc1 = converge(scene, &cfg, &r1);
    check("converges", rc1 == ALLOC_OK && r1.converged &&
                       r1.rmse <= cfg.target_rmse && r1.rounds <= MAX_ROUNDS);

    /* reproducibility: same seed + config must give bit-identical reruns */
    converge(scene, &cfg, &r2);
    check("converge_reproducible", memcmp(&r1, &r2, sizeof(r1)) == 0);

    long na[N], nb[N];
    double fa[N], fb[N];
    allocate_neyman(N, SIGMA, B, NMIN, na, fa, &minreq);
    allocate_neyman(N, SIGMA, B, NMIN, nb, fb, &minreq);
    check("alloc_reproducible",
          memcmp(na, nb, sizeof(na)) == 0 && memcmp(fa, fb, sizeof(fa)) == 0);

    /* variance estimates recover ground truth once every pixel is well sampled */
    allocate_neyman(N, SIGMA, (long)N * 256, NMIN, n, frac, &minreq);
    WStat st[N];
    Rng r;
    rng_seed(&r, 999u);
    memset(st, 0, sizeof(st));
    for (int i = 0; i < N; i++)
        for (long k = 0; k < n[i]; k++)
            stat_push(&st[i], pixel_sample(&scene[i], &r));
    int var_ok = 1;
    for (int i = 0; i < N; i++) {
        double rel = fabs(stat_var(&st[i]) - SIGMA[i] * SIGMA[i]) / (SIGMA[i] * SIGMA[i]);
        if (rel > 0.5) var_ok = 0;
    }
    check("variance_recovers", var_ok);

    /* empirical Monte Carlo: adaptive allocation beats uniform at equal budget */
    allocate_neyman(N, SIGMA, B, NMIN, n, frac, &minreq);
    double mse_a, mse_u;
    empirical_mse(scene, N, n, B / N, 400, 777u, &mse_a, &mse_u);
    double emp_ratio = sqrt(mse_a / mse_u);
    check("empirical_adaptive_wins", mse_a < mse_u * 0.9);
    check("empirical_matches_prediction", fabs(emp_ratio - err_ratio) < 0.15);

    /* ---------------- demo printout ---------------- */
    allocate_neyman(N, SIGMA, B, NMIN, n, frac, &minreq);
    printf("\n--- allocation demo: budget=%ld n_min=%d min_usable=%ld ---\n", B, NMIN, minreq);
    printf("pixel  sigma  samples  fraction\n");
    for (int i = 0; i < N; i++)
        printf("%5d  %5.2f  %7ld  %.4f\n", i, SIGMA[i], n[i], frac[i]);
    printf("rmse_uniform=%.6f rmse_adaptive=%.6f error_ratio=%.6f\n",
           predict_uniform_rmse(N, SIGMA, B), predict_rmse(N, SIGMA, n),
           predict_rmse(N, SIGMA, n) / predict_uniform_rmse(N, SIGMA, B));
    printf("converged=%d rounds=%d final_budget=%ld final_rmse=%.6f target=%.3f\n",
           r1.converged, r1.rounds, r1.budget, r1.rmse, cfg.target_rmse);
    int zrc = allocate_neyman(N, SIGMA, 0, NMIN, n, frac, &minreq);
    printf("zero_budget: rc=%d min_required=%ld\n", zrc, minreq);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
