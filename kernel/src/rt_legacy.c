/* rt_legacy.c -- the pre-unification random sources, one per site.
 *
 * This file exists only so the "before" behaviour can be measured and frozen
 * on this machine.  It is selected by building with -DRT_KERNEL_LEGACY=1.
 * Each site has its own generator, its own state and its own draw order, which
 * is exactly what made results machine- and schedule-dependent once work was
 * spread across threads.  Do not use this path for rendering.
 */
#include <stdint.h>
#include <math.h>
#include "rt_kernel.h"
#include "rt_shim.h"

/* Attenuation site's private LCG (numerical-recipes constants). */
static uint32_t s_atten_state = 0x12345678u;

/* Material site's private LCG (different constants, different state). */
static uint32_t s_material_state = 0x0BADF00Du;

/* Sampling site's private xorshift64. */
static uint64_t s_sampling_state = 0x9E3779B97F4A7C15ULL;

static uint32_t lcg_atten(void) {
    s_atten_state = s_atten_state * 1664525u + 1013904223u;
    return s_atten_state;
}
static uint32_t lcg_material(void) {
    s_material_state = s_material_state * 1103515245u + 12345u;
    return s_material_state;
}
static uint64_t xorshift_sampling(void) {
    uint64_t x = s_sampling_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    s_sampling_state = x;
    return x;
}

static double atten_u(void)  { return (double)(lcg_atten() >> 8) * (1.0 / 16777216.0); }
static float  material_u(void) { return (float)(lcg_material() >> 8) / 16777216.0f; }
static double sampling_u(void) { return (double)(xorshift_sampling() >> 11) * (1.0 / 9007199254740992.0); }

uint64_t rt_rng_configure(uint64_t legacy_seed) {
    s_atten_state = (uint32_t)(0x12345678u ^ legacy_seed);
    s_material_state = (uint32_t)(0x0BADF00Du ^ (legacy_seed >> 32));
    s_sampling_state = 0x9E3779B97F4A7C15ULL ^ legacy_seed;
    return legacy_seed;
}

RtAttenSample rt_atten_sample(double sigma, uint64_t key) {
    RtAttenSample r = { 0.0, 0.0 };
    (void)key; /* legacy path ignores the work-item key: shared state, ordered */
    double u = atten_u();
    if (u >= 1.0) u = 0.9999999999999999; /* keep -log1p finite */
    r.t = -log1p(-u) / sigma;
    r.weight = exp(-sigma * r.t);
    return r;
}

double rt_material_mean_cos(int n, uint64_t key) {
    (void)key;
    if (n <= 0) return 0.0;
    double acc = 0.0;
    float v[3];
    for (int i = 0; i < n; i++) {
        /* Sequenced explicitly: C leaves the evaluation order of function
           arguments unspecified, so two draws in one call would bind u1/u2
           differently under gcc and clang. */
        float u1 = material_u();
        float u2 = material_u();
        rt_shim_cosine_sample(u1, u2, v);
        acc += (double)v[2];
    }
    return acc / (double)n;
}

double rt_strat_quad(int strata, uint64_t key) {
    (void)key;
    if (strata <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < strata; i++) {
        double u = ((double)i + sampling_u()) / (double)strata;
        sum += rt_shim_f_quad(u);
    }
    return sum / (double)strata;
}
