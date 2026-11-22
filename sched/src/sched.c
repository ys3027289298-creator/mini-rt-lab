/* sched: tile-render load balancing by predicted cost, greedy placement,
   dynamic claiming fallback, and makespan/variance comparison. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_TILES 256
#define MAX_CORES 32

#define TILE_W 16
#define HOT_X 11.5
#define HOT_Y 3.2
#define HOT_SIGMA 2.2
#define HOT_GAIN 24.0

#define PILOT_SAMPLES 16
#define ERR_THRESHOLD 0.35

enum { MODE_STATIC = 0, MODE_GREEDY = 1, MODE_DYNAMIC = 2 };

typedef struct {
    int ntiles, ncores;
    double true_cost[MAX_TILES];
    double est_cost[MAX_TILES];
    double pilot_err; /* mean |est-true|/true on an evenly spaced pilot sample */
} Job;

typedef struct {
    int mode;
    int ntiles, ncores;
    int core_of[MAX_TILES];
    double finish[MAX_CORES];
    double makespan;
    double finish_var; /* population variance of per-core finish times */
    double utilization;
} Plan;

static double rnd01(unsigned* s) {
    *s = *s * 1664525u + 1013904223u;
    return (double)(*s >> 8) / 16777216.0;
}

static int job_valid(const Job* j) {
    return j && j->ntiles > 0 && j->ntiles <= MAX_TILES &&
           j->ncores > 0 && j->ncores <= MAX_CORES;
}

/* Synthetic 16-wide tile grid: one bright hotspot plus deterministic cost
   jitter. Estimator noise is multiplicative on the true cost. */
static int job_init(Job* j, int ntiles, int ncores, double noise, unsigned seed) {
    if (!j || ntiles <= 0 || ntiles > MAX_TILES ||
        ncores <= 0 || ncores > MAX_CORES || noise < 0.0) return -1;
    memset(j, 0, sizeof(*j));
    j->ntiles = ntiles;
    j->ncores = ncores;

    unsigned st = seed;
    for (int i = 0; i < ntiles; i++) {
        int x = i % TILE_W, y = i / TILE_W;
        double dx = (double)x - HOT_X, dy = (double)y - HOT_Y;
        double hot = HOT_GAIN * exp(-(dx * dx + dy * dy) / (2.0 * HOT_SIGMA * HOT_SIGMA));
        double base = 0.4 + 0.9 * rnd01(&st);
        j->true_cost[i] = base + hot;

        double f = 1.0 + noise * (2.0 * rnd01(&st) - 1.0);
        if (f < 0.05) f = 0.05;
        j->est_cost[i] = j->true_cost[i] * f;
    }

    int stride = ntiles / PILOT_SAMPLES;
    if (stride < 1) stride = 1;
    double acc = 0.0;
    int n = 0;
    for (int i = stride / 2; i < ntiles; i += stride) {
        acc += fabs(j->est_cost[i] - j->true_cost[i]) / j->true_cost[i];
        n++;
    }
    j->pilot_err = (n > 0) ? acc / n : 0.0;
    return 0;
}

static int plan_finalize(const Job* j, Plan* p, int mode) {
    p->mode = mode;
    p->ntiles = j->ntiles;
    p->ncores = j->ncores;
    double total = 0.0, max = 0.0;
    for (int c = 0; c < p->ncores; c++) {
        total += p->finish[c];
        if (p->finish[c] > max) max = p->finish[c];
    }
    p->makespan = max;
    double mean = total / p->ncores;
    double var = 0.0;
    for (int c = 0; c < p->ncores; c++) {
        double d = p->finish[c] - mean;
        var += d * d;
    }
    p->finish_var = var / p->ncores;
    p->utilization = (max > 0.0) ? total / ((double)p->ncores * max) : 0.0;
    return 0;
}

/* Static baseline: contiguous equal-count ranges, one per core. */
static int plan_static(const Job* j, Plan* p) {
    if (!job_valid(j) || !p) return -1;
    memset(p, 0, sizeof(*p));
    int base = j->ntiles / j->ncores;
    int rem = j->ntiles % j->ncores;
    int i = 0;
    for (int c = 0; c < j->ncores; c++) {
        int cnt = base + (c < rem ? 1 : 0);
        for (int k = 0; k < cnt; k++, i++) {
            p->core_of[i] = c;
            p->finish[c] += j->true_cost[i];
        }
    }
    return plan_finalize(j, p, MODE_STATIC);
}

/* LPT-style greedy: largest predicted tile first, onto the least predicted load. */
static int plan_greedy(const Job* j, Plan* p) {
    if (!job_valid(j) || !p) return -1;
    memset(p, 0, sizeof(*p));

    int order[MAX_TILES];
    for (int i = 0; i < j->ntiles; i++) order[i] = i;
    for (int i = 1; i < j->ntiles; i++) {
        int key = order[i], k = i - 1;
        while (k >= 0 && (j->est_cost[order[k]] < j->est_cost[key] ||
                          (j->est_cost[order[k]] == j->est_cost[key] && order[k] > key))) {
            order[k + 1] = order[k];
            k--;
        }
        order[k + 1] = key;
    }

    double est_load[MAX_CORES], true_load[MAX_CORES];
    memset(est_load, 0, sizeof(est_load));
    memset(true_load, 0, sizeof(true_load));
    for (int i = 0; i < j->ntiles; i++) {
        int idx = order[i], best = 0;
        for (int c = 1; c < j->ncores; c++)
            if (est_load[c] < est_load[best]) best = c;
        p->core_of[idx] = best;
        est_load[best] += j->est_cost[idx];
        true_load[best] += j->true_cost[idx];
    }
    for (int c = 0; c < j->ncores; c++) p->finish[c] = true_load[c];
    return plan_finalize(j, p, MODE_GREEDY);
}

/* Dynamic claiming: cores pull the next queued tile whenever they free up. */
static int plan_dynamic(const Job* j, Plan* p) {
    if (!job_valid(j) || !p) return -1;
    memset(p, 0, sizeof(*p));
    for (int i = 0; i < j->ntiles; i++) {
        int best = 0;
        for (int c = 1; c < j->ncores; c++)
            if (p->finish[c] < p->finish[best]) best = c;
        p->core_of[i] = best;
        p->finish[best] += j->true_cost[i];
    }
    return plan_finalize(j, p, MODE_DYNAMIC);
}

/* Auto mode: trust estimates only when the pilot sample says they are close. */
static int plan_auto(const Job* j, Plan* p) {
    if (!job_valid(j) || !p) return -1;
    if (j->pilot_err > ERR_THRESHOLD) return plan_dynamic(j, p);
    return plan_greedy(j, p);
}

static double plan_speedup(const Plan* baseline, const Plan* p) {
    if (!baseline || !p || p->makespan <= 0.0) return 0.0;
    return baseline->makespan / p->makespan;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int plan_covers_all(const Plan* p) {
    int seen[MAX_TILES];
    memset(seen, 0, sizeof(seen));
    for (int i = 0; i < p->ntiles; i++) {
        int c = p->core_of[i];
        if (c < 0 || c >= p->ncores) return 0;
        seen[i] = 1;
    }
    for (int i = 0; i < p->ntiles; i++) if (!seen[i]) return 0;
    return 1;
}

int main(void) {
    Job bad;
    Plan bp;
    check("zero_tiles_rejected", job_init(&bad, 0, 8, 0.05, 1) == -1);
    check("zero_cores_rejected", job_init(&bad, 64, 0, 0.05, 1) == -1);
    memset(&bad, 0, sizeof(bad));
    check("empty_plan_rejected", plan_auto(&bad, &bp) == -1);

    /* Accurate estimates: greedy placement should win big over static rows. */
    Job lo;
    check("job_low_ok", job_init(&lo, 256, 8, 0.05, 1234) == 0);
    Plan lo_static, lo_greedy, lo_dynamic, lo_auto;
    check("low_static_ok", plan_static(&lo, &lo_static) == 0);
    check("low_greedy_ok", plan_greedy(&lo, &lo_greedy) == 0);
    check("low_dynamic_ok", plan_dynamic(&lo, &lo_dynamic) == 0);
    check("low_auto_ok", plan_auto(&lo, &lo_auto) == 0);
    check("low_auto_picks_greedy", lo_auto.mode == MODE_GREEDY);
    check("low_covers_all", plan_covers_all(&lo_auto));

    double lo_speedup = plan_speedup(&lo_static, &lo_auto);
    printf("low_pilot_err=%.4f\n", lo.pilot_err);
    printf("low_static_makespan=%.3f\n", lo_static.makespan);
    printf("low_auto_makespan=%.3f\n", lo_auto.makespan);
    printf("low_speedup=%.3f\n", lo_speedup);
    printf("low_greedy_makespan=%.3f\n", lo_greedy.makespan);
    printf("low_dynamic_makespan=%.3f\n", lo_dynamic.makespan);
    printf("low_static_var=%.3f\n", lo_static.finish_var);
    printf("low_auto_var=%.3f\n", lo_auto.finish_var);
    check("low_beats_static", lo_auto.makespan < lo_static.makespan);
    check("low_speedup_gt_1_5", lo_speedup > 1.5);
    check("low_var_below_static", lo_auto.finish_var < lo_static.finish_var);
    check("low_util_above_0_9", lo_auto.utilization > 0.9);
    check("low_greedy_near_dynamic", lo_greedy.makespan <= lo_dynamic.makespan * 1.10);

    /* Noisy estimates: greedy misplaces hot tiles, so auto must claim dynamically. */
    Job hi;
    check("job_high_ok", job_init(&hi, 256, 8, 1.5, 1234) == 0);
    Plan hi_static, hi_greedy, hi_dynamic, hi_auto;
    check("high_static_ok", plan_static(&hi, &hi_static) == 0);
    check("high_greedy_ok", plan_greedy(&hi, &hi_greedy) == 0);
    check("high_dynamic_ok", plan_dynamic(&hi, &hi_dynamic) == 0);
    check("high_auto_ok", plan_auto(&hi, &hi_auto) == 0);
    check("high_auto_picks_dynamic", hi_auto.mode == MODE_DYNAMIC);
    check("high_covers_all", plan_covers_all(&hi_auto));

    double hi_speedup = plan_speedup(&hi_static, &hi_auto);
    printf("high_pilot_err=%.4f\n", hi.pilot_err);
    printf("high_static_makespan=%.3f\n", hi_static.makespan);
    printf("high_greedy_makespan=%.3f\n", hi_greedy.makespan);
    printf("high_dynamic_makespan=%.3f\n", hi_dynamic.makespan);
    printf("high_auto_makespan=%.3f\n", hi_auto.makespan);
    printf("high_speedup=%.3f\n", hi_speedup);
    check("high_beats_static", hi_auto.makespan < hi_static.makespan);
    check("high_speedup_gt_1_5", hi_speedup > 1.5);
    check("high_dynamic_beats_greedy", hi_dynamic.makespan < hi_greedy.makespan);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
