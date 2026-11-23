/* rt_unified.c -- all three sites draw from the one counter-based source.
 *
 * There is no mutable generator state anywhere in this file: every draw is a
 * pure function of (site, item, stream), so the value a work item receives is
 * independent of thread, scheduling order, process and machine.
 */
#include <math.h>
#include "rt_kernel.h"
#include "rt_rng.h"
#include "rt_shim.h"

uint64_t rt_rng_configure(uint64_t legacy_seed) {
    (void)legacy_seed; /* unified path has nothing to reseed */
    return 0;
}

RtAttenSample rt_atten_sample(double sigma, uint64_t key) {
    RtAttenSample r = { 0.0, 0.0 };
    double u = rt_rng_unit(key);
    r.t = -log1p(-u) / sigma;
    r.weight = exp(-sigma * r.t);
    return r;
}

double rt_material_mean_cos(int n, uint64_t key) {
    if (n <= 0) return 0.0;
    double acc = 0.0;
    float v[3];
    for (int i = 0; i < n; i++) {
        float u1 = (float)rt_rng_unit(rt_rng_key(RT_SITE_MATERIAL, key, (uint64_t)i * 2u));
        float u2 = (float)rt_rng_unit(rt_rng_key(RT_SITE_MATERIAL, key, (uint64_t)i * 2u + 1u));
        rt_shim_cosine_sample(u1, u2, v);
        acc += (double)v[2];
    }
    return acc / (double)n;
}

double rt_strat_quad(int strata, uint64_t key) {
    if (strata <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < strata; i++) {
        /* One stratified sample per stratum, jittered by the unified source:
           E[sum] = integral exactly, so the estimator stays unbiased. */
        double u = ((double)i + rt_rng_unit(rt_rng_key(RT_SITE_SAMPLING, key, (uint64_t)i)))
                   / (double)strata;
        sum += rt_shim_f_quad(u);
    }
    return sum / (double)strata;
}
