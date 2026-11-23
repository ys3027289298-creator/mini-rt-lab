/* scene_hash.c -- end-to-end fingerprint of one fixed scene.
 *
 * The scene is ITEMS work items; for every item the driver evaluates all
 * three sites (attenuation, material, sampling) and folds every raw output
 * double into a 64-bit FNV-1a hash.  One line is printed:
 *
 *   scene_hash <16 hex digits> items=<n> threads=<t> order=<seq|strided>
 *
 * Two builds (or two thread counts, or two orders) that produce the same
 * scene must print the same hash.  This is the "same scene on another
 * machine" check: a hash differing across runs or machines means some site's
 * random source is not a pure function of the work item.
 *
 * Usage: scene_hash [threads] [items] [order]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include "rt_kernel.h"

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

#define MAX_ITEMS 8192
#define MAX_THREADS 64

/* One scene item: all three sites, raw values kept for hashing. */
typedef struct { double atten, mat, strat; } Item;
static Item g_items_buf[MAX_ITEMS];
static int g_items;

static uint64_t fnv1a(uint64_t h, const void* data, size_t n) {
    const unsigned char* p = (const unsigned char*)data;
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}

static void eval_item(int i) {
    g_items_buf[i].atten = rt_atten_sample((double)PROBE_SIGMA, (uint64_t)i).weight;
    g_items_buf[i].mat   = rt_material_mean_cos(PROBE_N, (uint64_t)i);
    g_items_buf[i].strat = rt_strat_quad(PROBE_STRATA, (uint64_t)i);
}

typedef struct { int tid, threads; } Assign;
static void* worker(void* arg) {
    Assign* a = (Assign*)arg;
    rt_rng_configure(PROBE_SEED);
    for (int i = a->tid; i < g_items; i += a->threads) eval_item(i);
    return NULL;
}

int main(int argc, char** argv) {
    int threads = argc > 1 ? atoi(argv[1]) : 1;
    g_items     = argc > 2 ? atoi(argv[2]) : PROBE_ITEMS;
    const char* order = argc > 3 ? argv[3] : "seq";
    if (threads < 1) threads = 1;
    if (threads > MAX_THREADS) threads = MAX_THREADS;
    if (g_items < 0) g_items = 0;
    if (g_items > MAX_ITEMS) g_items = MAX_ITEMS;

    if (strcmp(order, "strided") == 0) {
        rt_rng_configure(PROBE_SEED);
        for (int t = 0; t < threads; t++)
            for (int i = t; i < g_items; i += threads) eval_item(i);
    } else if (threads == 1) {
        Assign one = { 0, 1 };
        worker(&one);
    } else {
        pthread_t th[MAX_THREADS];
        Assign as[MAX_THREADS];
        for (int t = 0; t < threads; t++) {
            as[t].tid = t; as[t].threads = threads;
            if (pthread_create(&th[t], NULL, worker, &as[t]) != 0) return 2;
        }
        for (int t = 0; t < threads; t++) pthread_join(th[t], NULL);
    }

    uint64_t h = 1469598103934665603ULL;
    h = fnv1a(h, &g_items, sizeof(g_items));
    for (int i = 0; i < g_items; i++) {
        h = fnv1a(h, &g_items_buf[i].atten, sizeof(double));
        h = fnv1a(h, &g_items_buf[i].mat,   sizeof(double));
        h = fnv1a(h, &g_items_buf[i].strat, sizeof(double));
    }
    printf("scene_hash %016llx items=%d threads=%d order=%s\n",
           (unsigned long long)h, g_items, threads, order);
    return 0;
}
