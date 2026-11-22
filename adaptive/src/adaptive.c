/* adaptive: per-pixel variance-driven early stopping for Monte Carlo integration.

   Two stopping rules share one sample stream (paired comparison):
     eb    : empirical Bernstein half-width, a *provable* high-probability
             bound with a uniform boundary over a fixed checkpoint schedule;
     naive : sample-SD normal approximation (no range term), the common
             "loose" heuristic that can under-cover small/heavy-tailed samples.

   Every pixel always draws n_max samples, so the full-sample reference is
   available for the same stream and early stopping never changes the RNG
   progression.  Ground-truth integrals are analytic, which lets every claimed
   bound be checked directly instead of asserted.
*/
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N_MAX        2048
#define N_PIXELS     2000
#define N_CRIT       2
#define N_SCEN       3
#define MIN_SAMPLES  8
#define GROW         1.25
#define TOL          0.03
#define ALPHA        0.05

/* ---------- PRNG: deterministic per-pixel streams ---------- */
static unsigned long long rng_state;
static void rng_seed(unsigned long long s) { rng_state = s ? s : 0x9e3779b97f4a7c15ULL; }
static double rng_next(void) {
    unsigned long long x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return (double)((x * 0x2545F4914F6CDD1DULL) >> 11) / (double)(1ULL << 53);
}

/* ---------- Welford online accumulator (unbiased variance) ---------- */
typedef struct { int n; double mean, M2; } Acc;
static void acc_push(Acc* a, double x) {
    a->n++;
    double d = x - a->mean;
    a->mean += d / a->n;
    a->M2 += d * (x - a->mean);
}
static double acc_var(const Acc* a) { return a->n > 1 ? a->M2 / (a->n - 1) : 0.0; }

/* ---------- checkpoint schedule (shared, fixed in advance) ---------- */
static int make_checkpoints(int* cp, int nmax) {
    int k = 0, prev = 0;
    for (double n = MIN_SAMPLES; n < nmax; n *= GROW) {
        int m = (int)n;
        if (m > prev) { cp[k++] = m; prev = m; }
    }
    cp[k++] = nmax;
    return k;
}

/* ---------- stopping rules ----------
   Empirical Bernstein (Maurer & Pontil 2009, unbiased-variance form):
       Pr( |mean - mu| > h ) <= 2 exp(- n h^2 / (2 s^2 + 2 R h/3))
   With range R and a union bound over k checkpoints, the half-width that
   keeps the per-pixel failure probability at most alpha is
       h = sqrt(2 s^2 log(2k/alpha)/n) + (2R/3) log(2k/alpha)/n.
*/
typedef struct {
    const char* name;
    double tol, alpha, R;
    int    eb;                 /* 1: empirical Bernstein, 0: naive Gaussian */
} Criterion;

static double halfwidth(const Criterion* c, double var, int n, int k) {
    if (c->eb) {
        double l = log(2.0 * k / c->alpha);
        return sqrt(2.0 * var * l / n) + (2.0 * c->R / 3.0) * l / n;
    }
    double z = 1.959964;        /* 0.975 quantile: nominal 95% only */
    return z * sqrt(var / n);
}

/* ---------- test integrands on [0,1): analytic truth and range ---------- */
static double f_quad(double u)  { return u * u; }
static double f_sine(double u)  { return sin(3.14159265358979 * u); }
static double f_spike(double u) { return u < 0.05 ? 1.0 : 0.0; }

typedef struct { const char* name; double (*f)(double); double truth, R; } Scene;

/* ---------- per-criterion aggregate metrics ---------- */
typedef struct {
    long   early, exhausted;
    long   early_viol;        /* early-stopped pixels with |est-truth| > tol */
    long   exhausted_viol;    /* same breach, but for pixels run to n_max */
    long   union_pixels;      /* pixels whose bound failed at >=1 checkpoint */
    long   diff_viol;         /* |early_est - full_est| > tol */
    double sum_n, max_n, min_n;
    double sum_abs_bias, max_abs_bias, sum_sq_bias;
    double sum_early_bias;    /* sum of |bias| over early-stopped pixels */
    double sum_cap_exh;       /* sum of n_max half-widths over exhausted pixels */
    double sum_full_sq;       /* full-sample MSE against truth */
    double sum_abs_diff, max_abs_diff, sum_sq_diff; /* vs full-sample estimate */
} Metrics;

static void m_reset(Metrics* m) { memset(m, 0, sizeof(*m)); m->min_n = N_MAX + 1; }
static void m_add(Metrics* m, int is_early, int stop_n, double est, double full,
                  double truth, double tol, double h_nmax) {
    double bias = fabs(est - truth), diff = fabs(est - full);
    m->sum_n += stop_n;
    if (stop_n > m->max_n) m->max_n = stop_n;
    if (stop_n < m->min_n) m->min_n = stop_n;
    m->sum_abs_bias += bias; if (bias > m->max_abs_bias) m->max_abs_bias = bias; m->sum_sq_bias += bias * bias;
    m->sum_abs_diff += diff; if (diff > m->max_abs_diff) m->max_abs_diff = diff; m->sum_sq_diff += diff * diff;
    m->sum_full_sq += (full - truth) * (full - truth);
    if (is_early) { m->early++; m->sum_early_bias += bias; if (bias > tol) m->early_viol++; }
    else          { m->exhausted++; m->sum_cap_exh += h_nmax; if (bias > tol) m->exhausted_viol++; }
    if (diff > tol) m->diff_viol++;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Run one pixel through both criteria on one shared, fully-consumed stream.
   Deterministic in (scene, pixel): the seed maps to a fixed xorshift stream.
   union_hit[c] is set when criterion c's half-width fails at any checkpoint
   for this pixel (pixel-level union-bound event). h_nmax[c] is the half-width
   evaluated at n_max, the good-event error cap for exhausted pixels. */
static void run_pixel(const Scene* sc, const Criterion* crit, int k, const int* cp,
                      int si, int p, int stop_n[N_CRIT], double est[N_CRIT],
                      int stopped[N_CRIT], int union_hit[N_CRIT],
                      double h_nmax[N_CRIT], double* full_out) {
    rng_seed(0xABCDEF01ULL + (unsigned long long)si * 1000003ULL + (unsigned long long)p * 97ULL);
    Acc acc = {0, 0.0, 0.0};
    for (int c = 0; c < N_CRIT; c++) { stopped[c] = 0; union_hit[c] = 0; }
    int cpi = 0;

    for (int n = 1; n <= N_MAX; n++) {
        acc_push(&acc, sc->f(rng_next()));
        if (cpi < k && n == cp[cpi]) {
            for (int c = 0; c < N_CRIT; c++) {
                if (stopped[c]) continue;
                double h = halfwidth(&crit[c], acc_var(&acc), n, k);
                if (fabs(acc.mean - sc->truth) > h) union_hit[c] = 1;
                if (n == N_MAX) h_nmax[c] = h;
                if (h <= crit[c].tol) {
                    stopped[c] = 1; stop_n[c] = n; est[c] = acc.mean;
                }
            }
            cpi++;
        }
    }
    *full_out = acc.mean;
    for (int c = 0; c < N_CRIT; c++)
        if (!stopped[c]) { stop_n[c] = N_MAX; est[c] = *full_out; }
}

int main(int argc, char** argv) {
    const char* csv_path = (argc > 1) ? argv[1] : "artifacts/adaptive_pixels.csv";
    FILE* csv = fopen(csv_path, "w");
    if (!csv) printf("csv_open_failed=%s\n", csv_path);
    else     fprintf(csv, "scene,criterion,pixel,early,stop_n,est,full,truth,bias,diff,viol_truth,viol_full\n");

    Scene scenes[N_SCEN] = {
        { "quad",  f_quad,  1.0 / 3.0,              1.0 },
        { "sine",  f_sine,  2.0 / 3.14159265358979, 1.0 },
        { "spike", f_spike, 0.05,                   1.0 },
    };
    Criterion crit[N_CRIT] = {
        { "eb",    TOL, ALPHA, 1.0, 1 },  /* provable 95% half-width */
        { "naive", TOL, ALPHA, 1.0, 0 },  /* 1.96 * SE, nominal only  */
    };

    int cp[64];
    int k = make_checkpoints(cp, N_MAX);
    printf("checkpoints=%d nmax=%d pixels=%d tol=%.3f alpha=%.2f\n", k, N_MAX, N_PIXELS, TOL, ALPHA);

    for (int si = 0; si < N_SCEN; si++) {
        Metrics m[N_CRIT];
        for (int c = 0; c < N_CRIT; c++) m_reset(&m[c]);

        for (int p = 0; p < N_PIXELS; p++) {
            int stop_n[N_CRIT], stopped[N_CRIT], union_hit[N_CRIT];
            double est[N_CRIT], full, h_nmax[N_CRIT];
            run_pixel(&scenes[si], crit, k, cp, si, p, stop_n, est, stopped,
                      union_hit, h_nmax, &full);
            for (int c = 0; c < N_CRIT; c++) {
                m_add(&m[c], stopped[c], stop_n[c], est[c], full,
                      scenes[si].truth, crit[c].tol, h_nmax[c]);
                m[c].union_pixels += union_hit[c];
                if (csv) fprintf(csv, "%s,%s,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d\n",
                    scenes[si].name, crit[c].name, p, stopped[c], stop_n[c], est[c], full,
                    scenes[si].truth, fabs(est[c] - scenes[si].truth), fabs(est[c] - full),
                    fabs(est[c] - scenes[si].truth) > crit[c].tol,
                    fabs(est[c] - full) > crit[c].tol);
            }
        }
        if (csv) fflush(csv);

        double full_rmse = sqrt(m[0].sum_full_sq / N_PIXELS);
        for (int c = 0; c < N_CRIT; c++) {
            double np = N_PIXELS;
            double mean_n = m[c].sum_n / np;
            double rmse = sqrt(m[c].sum_sq_bias / np);
            double rmse_diff = sqrt(m[c].sum_sq_diff / np);
            double early_frac = m[c].early / np;
            double early_viol_rate = m[c].early ? (double)m[c].early_viol / m[c].early : 0.0;
            printf("scene=%s crit=%s early=%.4f exhausted=%.4f stop_n_mean=%.1f stop_n_range=[%.0f,%.0f] "
                   "saves=%.4f bias_mean=%.6f bias_max=%.6f rmse=%.6f full_rmse=%.6f rmse_ratio=%.2f "
                   "diff_mean=%.6f diff_max=%.6f rmse_diff=%.6f "
                   "early_viol=%d(rate_among_early=%.4f) exhausted_viol=%d viol_full=%d(%.4f) union_pixels=%d(%.4f)\n",
                   scenes[si].name, crit[c].name,
                   early_frac, m[c].exhausted / np, mean_n, m[c].min_n, m[c].max_n,
                   1.0 - mean_n / N_MAX,
                   m[c].sum_abs_bias / np, m[c].max_abs_bias, rmse, full_rmse, rmse / full_rmse,
                   m[c].sum_abs_diff / np, m[c].max_abs_diff, rmse_diff,
                   (int)m[c].early_viol, early_viol_rate, (int)m[c].exhausted_viol,
                   (int)m[c].diff_viol, m[c].diff_viol / np,
                   (int)m[c].union_pixels, m[c].union_pixels / np);
        }

        /* ---- bound assertions, evaluated per scene ---- */
        double slack = 0.03;  /* binomial sampling slack over 2000 pixels */
        /* 1. Union bound: fraction of pixels with any checkpoint breach <= alpha. */
        check("eb_union_bound_held",
              (double)m[0].union_pixels / N_PIXELS <= ALPHA + slack);
        /* 2. Early-stop guarantee: among pixels the rule actually stopped,
              the fraction whose true error exceeds tol is <= alpha. */
        if (m[0].early > 0)
            check("eb_early_viol_within_alpha",
                  (double)m[0].early_viol / m[0].early <= ALPHA + slack);
        else
            printf("eb_early_viol_within_alpha=skip(no early stops)\n");
        /* 3a. Conditional expected bias of early stops:
              E|bias| | stop] <= tol + R*alpha (good event caps at tol,
              bad event costs at most R). */
        if (m[0].early > 0)
            check("eb_early_expected_bias_bound",
                  m[0].sum_early_bias / m[0].early <= TOL + crit[0].R * ALPHA + 1e-3);
        /* 3b. Same bound for exhausted pixels with tol replaced by the
              half-width actually reached at n_max. */
        if (m[0].exhausted > 0) {
            double cap = m[0].sum_cap_exh / m[0].exhausted + crit[0].R * ALPHA;
            double exh_bias = (m[0].sum_abs_bias - m[0].sum_early_bias) / m[0].exhausted;
            check("eb_exhausted_expected_bias_bound", exh_bias <= cap + 1e-3);
        }
        /* 4. EB must actually save samples somewhere, or the exercise is void. */
        check("eb_achieves_early_stops", m[0].early > 0);
        /* 5. The loose naive rule: any tol breach must be surfaced, not hidden.
              On spike it systematically under-covers; require it to be reported. */
        if (strcmp(scenes[si].name, "spike") == 0)
            check("naive_loose_threshold_reported",
                  m[1].early_viol + m[1].exhausted_viol > 0 &&
                  (double)m[1].early_viol / (m[1].early ? m[1].early : 1) > ALPHA);
        /* 6. EB's RMSE stays within 3x of the full-sample reference even
              though it uses fewer samples.  No such guarantee exists for naive;
              its rmse_ratio is printed above as evidence, not asserted. */
        check("eb_rmse_not_blown_up",
              sqrt(m[0].sum_sq_bias / N_PIXELS) < 3.0 * full_rmse + 1e-6);
    }
    if (csv) fclose(csv);

    /* ---- structural self-tests ---- */
    int cp2[64]; int k2 = make_checkpoints(cp2, N_MAX);
    int ordered = 1, grows = 1;
    for (int i = 1; i < k2; i++) {
        if (cp2[i] <= cp2[i-1]) ordered = 0;
        if (cp2[i-1] < MIN_SAMPLES) grows = 0;
    }
    check("checkpoints_ordered", ordered && grows && cp2[0] == MIN_SAMPLES && cp2[k2-1] == N_MAX);

    Acc t = {0, 0.0, 0.0};
    double data[5] = {2.0, 4.0, 4.0, 4.0, 6.0};   /* mean 4, unbiased var 2 */
    for (int i = 0; i < 5; i++) acc_push(&t, data[i]);
    check("welford_mean", fabs(t.mean - 4.0) < 1e-12);
    check("welford_var", fabs(acc_var(&t) - 2.0) < 1e-12);

    /* reproducibility: rerun scene 0 pixel 0 and require bit-identical output */
    int s1[N_CRIT], st1[N_CRIT], s2[N_CRIT], st2[N_CRIT];
    int uh1[N_CRIT], uh2[N_CRIT];
    double e1[N_CRIT], f1, e2[N_CRIT], f2, h1[N_CRIT], h2[N_CRIT];
    run_pixel(&scenes[0], crit, k2, cp2, 0, 0, s1, e1, st1, uh1, h1, &f1);
    run_pixel(&scenes[0], crit, k2, cp2, 0, 0, s2, e2, st2, uh2, h2, &f2);
    int repro = (f1 == f2);
    for (int c = 0; c < N_CRIT; c++)
        repro = repro && s1[c] == s2[c] && e1[c] == e2[c] && st1[c] == st2[c] && uh1[c] == uh2[c];
    check("rerun_bit_identical", repro);

    printf("csv=%s\n", csv_path);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
