
/* reflect: tiled reflection-map update with temporal accumulation.

   The reflection map is re-ray-traced in tiles. Each frame only K of T
   tiles are recomputed ("only half of the map is finished before submit"
   in the baseline policy, K = T/2); stale tiles keep last frame's
   result, and a refreshed tile blends the new trace with the old value
   (EMA temporal accumulation). A metallic object rotates at angular
   speed omega against a static environment, so the stale texels show up
   as reflection smear. Two independent metrics are measured:

     lag   - mean effective age (frames) of the displayed texels,
             tracked per texel with an age counter; the steady state has
             the analytic form lag = (1-alpha)I/alpha + (I-1)/2 with
             refresh interval I = T/K.
     smear - spatial offset of each tile vs the sharp ground truth,
             fitted per tile by SSE over candidate shifts (sub-texel
             parabolic refinement), reported in texels and radians.
             Expected to satisfy smear ~ omega * lag. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N 256                       /* reflection map texels            */
#define T 16                        /* tiles                            */
#define TILE (N / T)                /* texels per tile                  */
#define PI2 6.283185307179586

#define BASE_K 8                    /* baseline: half the map per frame */
#define ALPHA 0.5                   /* EMA blend factor (fixed)         */
#define LAG_LIMIT 2.0               /* frames; fast tier must beat this */
#define SMEAR_BUDGET 0.12           /* tolerated angular smear (rad)    */
#define SHIFT_MAX 40                /* per-tile SSE search (texels)     */
#define WARMUP 128
#define MEASURE 128

/* static environment seen in the reflection: wide smooth features so a
   16-texel tile still carries a locatable signal. */
static double env(double theta) {
    return 0.5 + 0.5 * cos(theta)
             + 0.35 * cos(2.0 * theta + 1.7)
             + 0.15 * cos(3.0 * theta + 0.4);
}

/* ground truth of texel i at frame f: object rotated by omega * f */
static double ground_truth(int i, int f, double omega) {
    double theta = PI2 * (double)i / (double)N - omega * (double)f;
    return env(theta);
}

/* predicted steady-state mean age for refresh every I frames */
static double model_lag(int interval) {
    double age0 = (1.0 - ALPHA) * (double)interval / ALPHA;
    return age0 + ((double)interval - 1.0) / 2.0;
}

typedef struct {
    double lag;          /* frames, from per-texel age counters   */
    double smear_tex;    /* texels, per-tile SSE image fit        */
    double smear_rad;    /* radians                               */
    double predicted;    /* analytic model lag, cross-check       */
    int updated[T];      /* tiles touched over the whole run      */
} Result;

static Result run(double omega, int K, int frames) {
    double map[N], age[N];
    Result r;
    memset(&r, 0, sizeof(r));
    for (int i = 0; i < N; i++) {
        map[i] = ground_truth(i, 0, omega);
        age[i] = 0.0;
    }

    double age_sum = 0.0, shift_sum = 0.0;
    long age_count = 0, shift_count = 0;

    for (int f = 1; f <= frames; f++) {
        int start = (f * K) % T;                     /* round-robin tiles */
        for (int tx = 0; tx < T; tx++) {
            int tile = (start + tx) % T;
            int do_update = tx < K;
            if (do_update && f > WARMUP) r.updated[tile]++;
            for (int j = 0; j < TILE; j++) {
                int i = tile * TILE + j;
                if (do_update) {
                    double x = ground_truth(i, f, omega);
                    map[i] = ALPHA * x + (1.0 - ALPHA) * map[i];
                    age[i] = (1.0 - ALPHA) * (age[i] + 1.0);
                } else {
                    age[i] += 1.0;                   /* reuse last frame  */
                }
            }
        }

        if (f > WARMUP) {
            for (int i = 0; i < N; i++) {
                age_sum += age[i];
                age_count++;
            }

            /* per-tile SSE: which shift of the sharp current truth best
               matches this (uniformly aged) tile of the displayed map */
            for (int tile = 0; tile < T; tile++) {
                double err[2 * SHIFT_MAX + 1];
                int base_i = tile * TILE;
                for (int s = -SHIFT_MAX; s <= SHIFT_MAX; s++) {
                    double e = 0.0;
                    for (int j = 0; j < TILE; j++) {
                        int i = base_i + j;
                        double d = map[i] - ground_truth((i + s + N) % N, f, omega);
                        e += d * d;
                    }
                    err[s + SHIFT_MAX] = e;
                }
                int best = 0;
                for (int k = 1; k <= 2 * SHIFT_MAX; k++)
                    if (err[k] < err[best]) best = k;
                double sub = 0.0;
                if (best > 0 && best < 2 * SHIFT_MAX) {
                    double em = err[best - 1], e0 = err[best], ep = err[best + 1];
                    double den = em - 2.0 * e0 + ep;
                    if (fabs(den) > 1e-12) sub = 0.5 * (em - ep) / den;
                }
                shift_sum += (double)(best - SHIFT_MAX) + sub;
                shift_count++;
            }
        }
    }

    r.lag = age_sum / (double)age_count;
    r.smear_tex = shift_sum / (double)shift_count;
    r.smear_rad = r.smear_tex * PI2 / (double)N;
    r.predicted = model_lag(T / K);
    return r;
}

/* adaptive policy: pick the smallest per-frame tile count whose model
   lag keeps the angular smear omega * lag within the budget */
static int adapt_K(double omega) {
    static const int ks[5] = { 1, 2, 4, 8, 16 };
    for (int j = 0; j < 5; j++) {
        if (omega * model_lag(T / ks[j]) <= SMEAR_BUDGET) return ks[j];
    }
    return 16;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    static const double omegas[3] = { 0.005, 0.03, 0.06 };
    static const char* tier[3] = { "slow", "medium", "fast" };
    int total = WARMUP + MEASURE;

    /* stationary object: temporal accumulation must not create smear */
    Result st = run(0.0, BASE_K, total);
    printf("== stationary object (omega=0) ==\n");
    printf("static_lag_frames=%.3f static_smear_tex=%.3f static_smear_rad=%.4f\n",
           st.lag, st.smear_tex, st.smear_rad);

    printf("\n== baseline: fixed split r=1/2 (K=%d/%d), alpha=%.2f ==\n",
           BASE_K, T, ALPHA);
    printf("tier   omega   ratio  lag_frames  smear_tex  smear_rad\n");
    Result base[3];
    for (int k = 0; k < 3; k++) {
        base[k] = run(omegas[k], BASE_K, total);
        printf("%-6s %-6.3f 1/2    %9.3f  %9.3f  %9.4f\n",
               tier[k], omegas[k], base[k].lag,
               base[k].smear_tex, base[k].smear_rad);
    }

    int triggered = base[2].lag > LAG_LIMIT;
    printf("\nfast lag %.3f frames vs limit %.1f -> fixed_ratio_triggered=%d\n",
           base[2].lag, LAG_LIMIT, triggered);

    printf("\n== adaptive: split ratio vs lag (smear budget %.2f rad) ==\n",
           SMEAR_BUDGET);
    printf("tier   omega   ratio  K/T  alpha  lag_frames  smear_tex  smear_rad\n");
    Result adapt[3];
    int Ka[3];
    for (int k = 0; k < 3; k++) {
        Ka[k] = adapt_K(omegas[k]);
        adapt[k] = run(omegas[k], Ka[k], total);
        printf("%-6s %-6.3f %2d/%-2d  %.2f   %9.3f  %9.3f  %9.4f\n",
               tier[k], omegas[k], Ka[k], T, ALPHA, adapt[k].lag,
               adapt[k].smear_tex, adapt[k].smear_rad);
    }

    printf("\nchecks:\n");
    check("static_no_smear", fabs(st.smear_tex) < 0.01);
    int cover = 1;
    for (int tx = 0; tx < T; tx++)
        if (st.updated[tx] != MEASURE * BASE_K / T) cover = 0;
    check("tile_coverage_round_robin", cover);

    int model_ok = 1, smear_model_ok = 1, mono = 1;
    for (int k = 0; k < 3; k++) {
        if (fabs(base[k].lag - base[k].predicted) > 0.3) model_ok = 0;
        if (fabs(adapt[k].lag - adapt[k].predicted) > 0.3) model_ok = 0;
        double expect_tex = omegas[k] * base[k].lag * N / PI2;
        double tol = 0.15 * expect_tex + 0.2;   /* EMA blur confounds fast fits */
        if (fabs(base[k].smear_tex - expect_tex) > tol) smear_model_ok = 0;
        if (k > 0 && base[k].smear_rad <= base[k - 1].smear_rad) mono = 0;
    }
    check("lag_matches_age_model", model_ok);
    check("smear_matches_omega_lag", smear_model_ok);
    check("smear_grows_with_speed", mono);
    check("fast_baseline_over_limit", triggered);
    check("fast_adaptive_under_limit", adapt[2].lag < LAG_LIMIT);

    int budget_ok = 1;
    for (int k = 0; k < 3; k++)
        if (omegas[k] * adapt[k].lag > SMEAR_BUDGET + 1e-9) budget_ok = 0;
    check("adaptive_within_smear_budget", budget_ok);
    check("ratios_follow_angular_rate", Ka[0] < Ka[1] && Ka[1] < Ka[2]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
