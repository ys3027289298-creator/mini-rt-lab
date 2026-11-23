/* rt_kernel.h -- frozen public API of the rendering kernel.
 *
 * These declarations are the contract: signatures must not change when the
 * random source is unified.  `kernel/scripts/check.sh` verifies this header
 * against baseline/public_api.sha256 as well as with compile-time function
 * pointer assertions (tools/sigcheck.c) and an nm symbol diff between the
 * legacy and unified builds.
 */
#ifndef RT_KERNEL_H
#define RT_KERNEL_H

#include <stdint.h>

/* One free-path attenuation sample: distance t and its compensated weight. */
typedef struct {
    double t;
    double weight;
} RtAttenSample;

/* Attenuation site: sample t ~ Exp(sigma) and return exp(-sigma*t).
   Requires sigma > 0 and finite; on invalid input returns {0, 0}.
   E[weight] = 1/2 for every sigma > 0. */
RtAttenSample rt_atten_sample(double sigma, uint64_t key);

/* Material site: mean cosine of `n` cosine-weighted hemisphere samples,
   using the material module's cosine_sample().  E[mean cos] -> 2/3.
   n <= 0 returns 0. */
double rt_material_mean_cos(int n, uint64_t key);

/* Sampling site: jittered stratified estimate of integral_0^1 u^2 du,
   using the sampling module's f_quad().  Exactly unbiased; value 1/3.
   strata <= 0 returns 0. */
double rt_strat_quad(int strata, uint64_t key);

/* Additive configuration hook (not part of the frozen stochastic API):
   seeds the legacy per-site RNGs.  Has no effect in the unified build.
   Returns the seed actually used. */
uint64_t rt_rng_configure(uint64_t legacy_seed);

#endif /* RT_KERNEL_H */
