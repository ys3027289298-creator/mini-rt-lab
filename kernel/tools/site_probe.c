/* site_probe.c -- prints one site's output for a range of work items.
 *
 * Shared by the baseline capture, the legacy-vs-unified comparison, the
 * thread-independence check and the sanitizer runs.  The site and its
 * parameters are selected by macros so the same source can be built for every
 * combination:
 *
 *   -DPROBE_SITE=0  attenuation,      -DPROBE_SIGMA=...
 *   -DPROBE_SITE=1  material mean cos, -DPROBE_N=...
 *   -DPROBE_SITE=2  stratified quad,   -DPROBE_STRATA=...
 *
 * Usage: site_probe [threads] [items] [order]
 *        order = seq (default) | strided
 *
 * "seq" evaluates items 0..n-1 in order on one worker when threads == 1; with
 * threads > 1 items are strided across real threads (i -> i % threads), which
 * exercises cross-thread interleaving of a shared generator.
 *
 * "strided" is a single-threaded, race-free emulation of that interleaving: it
 * visits items in the order thread 0's share, then thread 1's share, ... and
 * stores results by item index.  A schedule-dependent generator produces
 * different values under this order; a pure counter-based source does not.
 * This makes the defect reproducible without relying on an actual data race.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include "rt_kernel.h"

#ifndef PROBE_SITE
#define PROBE_SITE 0
#endif
#ifndef PROBE_ITEMS
#define PROBE_ITEMS 256
#endif
#ifndef PROBE_SIGMA
#define PROBE_SIGMA 1.7
#endif
#ifndef PROBE_N
#define PROBE_N 8
#endif
#ifndef PROBE_STRATA
#define PROBE_STRATA 16
#endif
#ifndef PROBE_SEED
#define PROBE_SEED 0x5EEDu
#endif

#define PROBE_MAX_ITEMS 8192
#define PROBE_MAX_THREADS 64

static double g_out[PROBE_MAX_ITEMS];
static int g_items;

static double eval_site(int item) {
    switch (PROBE_SITE) {
    case 0: return rt_atten_sample((double)PROBE_SIGMA, (uint64_t)item).weight;
    case 1: return rt_material_mean_cos(PROBE_N, (uint64_t)item);
    default: return rt_strat_quad(PROBE_STRATA, (uint64_t)item);
    }
}

typedef struct { int tid, threads; } Assign;

static void* worker(void* arg) {
    Assign* a = (Assign*)arg;
    rt_rng_configure(PROBE_SEED);
    for (int i = a->tid; i < g_items; i += a->threads) g_out[i] = eval_site(i);
    return NULL;
}

int main(int argc, char** argv) {
    int threads = argc > 1 ? atoi(argv[1]) : 1;
    g_items    = argc > 2 ? atoi(argv[2]) : PROBE_ITEMS;
    const char* order = argc > 3 ? argv[3] : "seq";
    if (threads < 1) threads = 1;
    if (threads > PROBE_MAX_THREADS) threads = PROBE_MAX_THREADS;
    if (g_items < 0) g_items = 0;
    if (g_items > PROBE_MAX_ITEMS) g_items = PROBE_MAX_ITEMS;

    if (strcmp(order, "strided") == 0) {
        /* Sequential emulation of the thread-strided visit order. */
        rt_rng_configure(PROBE_SEED);
        for (int t = 0; t < threads; t++)
            for (int i = t; i < g_items; i += threads) g_out[i] = eval_site(i);
    } else if (threads == 1) {
        Assign one = { 0, 1 };
        (void)worker(&one);
    } else {
        pthread_t th[PROBE_MAX_THREADS];
        Assign as[PROBE_MAX_THREADS];
        for (int t = 0; t < threads; t++) {
            as[t].tid = t;
            as[t].threads = threads;
            if (pthread_create(&th[t], NULL, worker, &as[t]) != 0) return 2;
        }
        for (int t = 0; t < threads; t++) pthread_join(th[t], NULL);
    }

    for (int i = 0; i < g_items; i++) printf("%d %.17g\n", i, g_out[i]);
    return 0;
}
