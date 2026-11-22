
/* adaptive: block-wise adaptive sampling.
 *
 * Scene model: the [0,1)^2 integration domain is split into BX*BY blocks.
 * A disk ("specular highlight") produces rare bright draws, so blocks inside
 * it have high variance and uniform sampling leaves clustered noise there.
 * Each draw is Bernoulli: V_HI with probability p(x,y), otherwise V_LO.
 *
 * Pipeline:
 *   1. pilot pass with online Welford mean/variance per block
 *   2. allocate the remaining budget proportional to block variance
 *   3. invert a target standard error: n = (z*sigma/eps)^2
 *      (eps <= 0 or NaN is rejected: no finite n guarantees zero error)
 *   4. run uniform vs adaptive at an equal total sample count, compare
 *      total / worst-block / highlight-region MSE and give the ratios
 *   5. emit an ASCII error-distribution map and PGM heatmaps; report every
 *      block whose standard error stays above the target
 *
 * Note on the metric: n proportional to variance (n ~ sigma^2) equalizes the
 * per-block squared error, which removes the noisy hot spots, but at a fixed
 * total budget it does not lower the *mean* MSE -- the MSE-optimal split is
 * Neyman allocation n ~ sigma. The report therefore shows total, worst-block
 * and highlight-region ratios separately instead of hiding this behind one
 * flattering number.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#define BX 8
#define BY 8
#define NB (BX * BY)
#define PILOT 32                /* pilot samples per block */
#define EXTRA 4096              /* budget distributed after the pilot */
#define TOTAL (PILOT * NB + EXTRA)
#define REPS 8                  /* independent repeats for the comparison */

#define DISK_X 0.5
#define DISK_Y 0.5
#define DISK_R 0.18
#define P_IN 0.30               /* bright-draw probability inside the disk */
#define P_OUT 0.01              /* bright-draw probability outside the disk */
#define V_HI 8.0
#define V_LO 0.1

/* ---------------- deterministic RNG: splitmix seeding + xorshift64* ------- */

static uint64_t splitmix64(uint64_t* x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static uint64_t block_stream(int seed0, int b) {
    uint64_t z = (uint64_t)seed0 + (uint64_t)b * 0x9e3779b97f4a7c15ULL;
    splitmix64(&z);
    return splitmix64(&z);
}

static uint64_t rng_next(uint64_t* s) {
    uint64_t x = *s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    *s = x;
    return x * 0x2545f4914f6cdd1dULL;
}

static double rng_unit(uint64_t* s) {
    return (double)(rng_next(s) >> 11) * (1.0 / 9007199254740992.0);
}

/* ---------------- scene ---------------------------------------------------- */

static void block_center(int b, double* x, double* y) {
    int bx = b % BX, by = b / BX;
    *x = ((double)bx + 0.5) / BX;
    *y = ((double)by + 0.5) / BY;
}

static int in_disk(double x, double y) {
    double dx = x - DISK_X, dy = y - DISK_Y;
    return dx * dx + dy * dy <= DISK_R * DISK_R;
}

static double scene_p(double x, double y) {
    return in_disk(x, y) ? P_IN : P_OUT;
}

static double scene_mean(double p) { return p * V_HI + (1.0 - p) * V_LO; }
static double scene_var(double p) {
    double d = V_HI - V_LO;
    return p * (1.0 - p) * d * d;
}

static double scene_draw(int b, uint64_t* s) {
    double x, y;
    block_center(b, &x, &y);
    return rng_unit(s) < scene_p(x, y) ? V_HI : V_LO;
}

/* ---------------- online block statistics (Welford) ----------------------- */

typedef struct { double mean, m2, var; int n; } BlockStat;

static void welford_push(BlockStat* st, double v) {
    st->n++;
    double d = v - st->mean;
    st->mean += d / st->n;
    st->m2 += d * (v - st->mean);
}

static void welford_finish(BlockStat* st) {
    st->var = (st->n > 1) ? st->m2 / (st->n - 1) : 0.0;
    if (st->var < 0.0) st->var = 0.0;
}

/* ---------------- target-error inversion ---------------------------------- */

/* Samples needed for standard error eps with confidence factor z:
 *   n = var / (eps/z)^2 = (z*sigma/eps)^2.
 * Returns -1 when the request is impossible / ill-posed (eps <= 0, NaN,
 * bad z), 0 when the block has no variance. */
static long long required_samples(double var, double eps, double z) {
    if (!(eps > 0.0) || !(z > 0.0) || isnan(var)) return -1;
    if (var <= 0.0) return 0;
    double e = eps / z;
    double need = var / (e * e);
    if (!(need <= 9.0e17)) return 1000000000000000000LL;
    long long n = (long long)ceil(need);
    return n < 1 ? 1 : n;
}

/* ---------------- variance-proportional allocation ------------------------ */

/* Distributes exactly `extra` samples with n_i proportional to var_i.
 * Integer parts are handed out first, leftovers go to the largest fractional
 * parts (Hare quota). With no signal anywhere, fall back to round robin so
 * every block still receives samples. */
static void allocate(int extra, const double* var, int* alloc) {
    double sum = 0.0;
    for (int i = 0; i < NB; i++) sum += var[i];

    if (!(sum > 0.0)) {
        for (int i = 0; i < NB; i++) alloc[i] = 0;
        for (int k = 0; k < extra; k++) alloc[k % NB]++;
        return;
    }

    double frac[NB];
    int given = 0;
    for (int i = 0; i < NB; i++) {
        double exact = (double)extra * var[i] / sum;
        alloc[i] = (int)floor(exact);
        frac[i] = exact - (double)alloc[i];
        given += alloc[i];
    }
    while (given < extra) {
        int best = -1;
        double best_frac = -1.0;
        for (int i = 0; i < NB; i++) {
            if (frac[i] > best_frac) { best_frac = frac[i]; best = i; }
        }
        alloc[best]++;
        frac[best] = -1.0;
        given++;
    }
}

/* ---------------- sampling runs ------------------------------------------- */

static void run_uniform(int per_block, int seed0, BlockStat* out) {
    for (int b = 0; b < NB; b++) {
        uint64_t s = block_stream(seed0, b);
        memset(&out[b], 0, sizeof(out[b]));
        for (int k = 0; k < per_block; k++) welford_push(&out[b], scene_draw(b, &s));
        welford_finish(&out[b]);
    }
}

static void run_adaptive(const int* ns, int seed0, BlockStat* out) {
    for (int b = 0; b < NB; b++) {
        uint64_t s = block_stream(seed0, b);
        memset(&out[b], 0, sizeof(out[b]));
        for (int k = 0; k < ns[b]; k++) welford_push(&out[b], scene_draw(b, &s));
        welford_finish(&out[b]);
    }
}

/* Per-block squared error against the analytic reference means, plus the
 * total / highlight / worst-case aggregates. */
typedef struct { double total, hot, worst; int nhot; } Metrics;

static Metrics measure(const BlockStat* st) {
    Metrics m;
    memset(&m, 0, sizeof(m));
    double hot_acc = 0.0;
    for (int b = 0; b < NB; b++) {
        double x, y;
        block_center(b, &x, &y);
        double ref = scene_mean(scene_p(x, y));
        double d = st[b].mean - ref;
        double se2 = d * d;
        m.total += se2;
        if (se2 > m.worst) m.worst = se2;
        if (in_disk(x, y)) { hot_acc += se2; m.nhot++; }
    }
    m.total /= NB;
    m.hot = m.nhot ? hot_acc / m.nhot : 0.0;
    return m;
}

/* ---------------- error-distribution maps --------------------------------- */

static int write_pgm(const char* path, const double* v) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;
    double mx = 1e-12;
    for (int i = 0; i < NB; i++) if (v[i] > mx) mx = v[i];
    fprintf(f, "P5\n%d %d\n255\n", BX, BY);
    for (int by = 0; by < BY; by++) {
        for (int bx = 0; bx < BX; bx++) {
            double t = sqrt(v[by * BX + bx] / mx);
            int g = (int)floor(255.0 * t);
            if (g < 0) g = 0;
            if (g > 255) g = 255;
            fputc((unsigned char)g, f);
        }
    }
    fclose(f);
    return 1;
}

static const char MAP_CHARS[] = " .:-=+*#%@";

static void ascii_map(const char* title, const double* v) {
    double mx = 1e-12;
    for (int i = 0; i < NB; i++) if (v[i] > mx) mx = v[i];
    printf("%s (dark=low, '@'=max %.4f)\n", title, mx);
    for (int by = BY - 1; by >= 0; by--) {
        printf("  |");
        for (int bx = 0; bx < BX; bx++) {
            int g = (int)(9.0 * sqrt(v[by * BX + bx] / mx));
            if (g > 9) g = 9;
            putchar(MAP_CHARS[g]);
        }
        printf("|\n");
    }
}

static int cmp_double(const void* a, const void* b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

static double quantile(double* sorted, int n, double q) {
    double pos = q * (n - 1);
    int lo = (int)pos;
    int hi = lo + 1 < n ? lo + 1 : lo;
    double t = pos - lo;
    return sorted[lo] * (1.0 - t) + sorted[hi] * t;
}

/* ---------------- test harness --------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ===== 1. target-error inversion, zero target must be rejected ===== */
    printf("== target error inversion ==\n");
    long long rz = required_samples(1.0, 0.0, 1.0);
    printf("eps=0.0 -> %lld (rejected)\n", rz);
    check("reject_zero_target", rz == -1);
    check("reject_negative_target", required_samples(1.0, -0.5, 1.0) == -1);
    check("reject_nan_target", required_samples(1.0, NAN, 1.0) == -1);
    check("reject_bad_z", required_samples(1.0, 0.1, 0.0) == -1);
    check("invert_basic", required_samples(1.0, 0.5, 1.0) == 4);
    check("invert_zfactor", required_samples(1.0, 0.5, 2.0) == 16);
    check("zero_variance_free", required_samples(0.0, 0.1, 1.0) == 0);

    /* ===== 2. Welford against a hand-computed sequence: 1..5 ========== */
    BlockStat w;
    memset(&w, 0, sizeof(w));
    for (int k = 1; k <= 5; k++) welford_push(&w, (double)k);
    welford_finish(&w);
    check("welford_mean", fabs(w.mean - 3.0) < 1e-12);
    check("welford_var", fabs(w.var - 2.5) < 1e-12);

    /* ===== 3. variance-proportional allocation on oracle variances ===== */
    printf("\n== variance-proportional allocation ==\n");
    double oracle_var[NB];
    int oracle_alloc[NB], hot = -1, cold = -1;
    for (int b = 0; b < NB; b++) {
        double x, y;
        block_center(b, &x, &y);
        oracle_var[b] = scene_var(scene_p(x, y));
        if (in_disk(x, y)) { if (hot < 0) hot = b; }
        else if (cold < 0) cold = b;
    }
    allocate(EXTRA, oracle_var, oracle_alloc);
    int sum_alloc = 0, min_alloc = EXTRA, max_alloc = 0;
    for (int b = 0; b < NB; b++) {
        sum_alloc += oracle_alloc[b];
        if (oracle_alloc[b] < min_alloc) min_alloc = oracle_alloc[b];
        if (oracle_alloc[b] > max_alloc) max_alloc = oracle_alloc[b];
    }
    check("alloc_exhausts_budget", sum_alloc == EXTRA);
    check("alloc_highlight_gets_more", oracle_alloc[hot] > oracle_alloc[cold]);
    printf("redistributed %d extra samples: per block min=%d max=%d "
           "(hot block %d vs cold block %d)\n",
           EXTRA, min_alloc, max_alloc, oracle_alloc[hot], oracle_alloc[cold]);

    double zero_var[NB];
    int zero_alloc[NB];
    memset(zero_var, 0, sizeof(zero_var));
    allocate(EXTRA, zero_var, zero_alloc);
    int zero_sum = 0, zero_nonempty = 1;
    for (int b = 0; b < NB; b++) {
        zero_sum += zero_alloc[b];
        if (zero_alloc[b] == 0) zero_nonempty = 0;
    }
    check("alloc_no_signal_roundrobin", zero_sum == EXTRA && zero_nonempty);

    /* Analytic prediction of what n~sigma^2 achieves at this total budget.
     * It equalizes per-block SE^2 to S/N, so total MSE is unchanged but the
     * hot / worst-block error collapses -- the actual goal here. */
    double S = 0.0, v_hot = 0.0;
    for (int b = 0; b < NB; b++) S += oracle_var[b];
    v_hot = oracle_var[hot];
    int n_uniform = TOTAL / NB;
    double eq_se2 = S / TOTAL;                                  /* adaptive */
    double u_hot_se2 = v_hot / n_uniform;                       /* uniform  */
    double pred_hot_ratio = eq_se2 / u_hot_se2;
    double pred_total_ratio = eq_se2 / (S / NB / n_uniform);   /* == 1.0 */
    printf("oracle prediction: equalized SE^2=%.4f, total-MSE ratio=%.3f, "
           "hot-block ratio=%.3f\n",
           eq_se2, pred_total_ratio, pred_hot_ratio);
    check("oracle_hot_block_improves", pred_hot_ratio < 0.3);
    check("oracle_total_is_redistribution", fabs(pred_total_ratio - 1.0) < 1e-9);

    /* ===== 4. equal-budget uniform vs adaptive, repeated =============== */
    printf("\n== equal-budget comparison (%d repeats, %d samples each) ==\n",
           REPS, TOTAL);
    double tu = 0, hu = 0, wu = 0, ta = 0, ha = 0, wa = 0;
    BlockStat prod_u[NB], prod_a[NB];
    int prod_ns[NB];
    for (int r = 0; r < REPS; r++) {
        int seed0 = 1000 + r * 7919;

        BlockStat pilot[NB];
        double pilot_var[NB];
        for (int b = 0; b < NB; b++) {
            uint64_t s = block_stream(seed0, b);
            memset(&pilot[b], 0, sizeof(pilot[b]));
            for (int k = 0; k < PILOT; k++) welford_push(&pilot[b], scene_draw(b, &s));
            welford_finish(&pilot[b]);
            pilot_var[b] = pilot[b].var;
        }
        int extra_alloc[NB], ns[NB], used = 0;
        allocate(EXTRA, pilot_var, extra_alloc);
        for (int b = 0; b < NB; b++) { ns[b] = PILOT + extra_alloc[b]; used += ns[b]; }
        if (r == 0) check("adaptive_matches_budget", used == TOTAL);

        BlockStat u[NB], a[NB];
        run_uniform(n_uniform, seed0, u);
        run_adaptive(ns, seed0, a);
        Metrics mu = measure(u), ma = measure(a);
        tu += mu.total; hu += mu.hot; wu += mu.worst;
        ta += ma.total; ha += ma.hot; wa += ma.worst;

        if (r == REPS - 1) {
            memcpy(prod_u, u, sizeof(u));
            memcpy(prod_a, a, sizeof(a));
            memcpy(prod_ns, ns, sizeof(ns));
        }
    }
    tu /= REPS; hu /= REPS; wu /= REPS;
    ta /= REPS; ha /= REPS; wa /= REPS;
    printf("                 uniform    adaptive   ratio(adaptive/uniform)\n");
    printf("  total MSE      %.5f   %.5f   %.3f\n", tu, ta, ta / tu);
    printf("  highlight MSE  %.5f   %.5f   %.3f\n", hu, ha, ha / hu);
    printf("  worst-block    %.5f   %.5f   %.3f\n", wu, wa, wa / wu);
    check("adaptive_cuts_highlight_noise", ha / hu < 0.7);
    check("adaptive_cuts_worst_block", wa / wu < 0.7);
    check("adaptive_keeps_total_in_band", ta / tu < 1.5);

    /* ===== 5. error-distribution maps of the last production run ======== */
    printf("\n== error distribution (last run) ==\n");
    double se_u[NB], se_a[NB], ns_d[NB];
    double se_sorted[NB];
    for (int b = 0; b < NB; b++) {
        double x, y;
        block_center(b, &x, &y);
        double ref = scene_mean(scene_p(x, y));
        se_u[b] = fabs(prod_u[b].mean - ref);
        se_a[b] = fabs(prod_a[b].mean - ref);
        se_sorted[b] = se_a[b];
        ns_d[b] = (double)prod_ns[b];
    }
    ascii_map("uniform |error|", se_u);
    ascii_map("adaptive |error|", se_a);
    ascii_map("adaptive samples/block", ns_d);
    qsort(se_sorted, NB, sizeof(double), cmp_double);
    printf("adaptive |error| quantiles: p50=%.4f p90=%.4f p99=%.4f max=%.4f\n",
           quantile(se_sorted, NB, 0.50), quantile(se_sorted, NB, 0.90),
           quantile(se_sorted, NB, 0.99), se_sorted[NB - 1]);
    int ok_se = write_pgm("error_se.pgm", se_a);
    int ok_ns = write_pgm("samples_n.pgm", ns_d);
    check("pgm_error_map_written", ok_se);
    check("pgm_sample_map_written", ok_ns);
    printf("heatmaps written: error_se.pgm, samples_n.pgm\n");

    /* ===== 6. target error: convergence report ========================== */
    printf("\n== target-error convergence (adaptive run) ==\n");
    double se_est[NB];
    for (int b = 0; b < NB; b++)
        se_est[b] = (prod_a[b].n > 1) ? sqrt(prod_a[b].var / prod_a[b].n) : 0.0;

    double targets[2] = { 0.10, 0.30 };
    for (int t = 0; t < 2; t++) {
        double eps = targets[t];
        int bad = 0;
        long long need_total = 0;
        for (int b = 0; b < NB; b++) {
            long long need = required_samples(oracle_var[b], eps, 1.0);
            if (need > 0) need_total += need;
            if (se_est[b] > eps) bad++;
        }
        printf("target eps=%.2f: unconverged blocks=%d/%d, total samples "
               "required for target=%lld (budget=%d)\n",
               eps, bad, NB, need_total, TOTAL);
        if (bad > 0) {
            printf("  first unconverged (block x,y: used_n, est_SE, required_n):\n");
            int shown = 0;
            for (int b = 0; b < NB && shown < 6; b++) {
                if (se_est[b] <= eps) continue;
                long long need = required_samples(oracle_var[b], eps, 1.0);
                printf("   (%d,%d): n=%d SE=%.3f need=%lld\n",
                       b % BX, b / BX, prod_a[b].n, se_est[b], need);
                shown++;
            }
        }
        if (t == 0) check("tight_target_reports_unconverged", bad > 0);
        else check("loose_target_all_converged", bad == 0);
    }

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
