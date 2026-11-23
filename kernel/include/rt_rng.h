/* rt_rng.h -- the single random source of the rendering kernel.
 *
 * Counter-based: every draw is a pure function of a 64-bit key.  There is no
 * mutable generator state, so the value assigned to a work item does not
 * depend on which thread runs it, in what order, or on which machine it is
 * built.  Different sites and samples are domain-separated through the key.
 */
#ifndef RT_RNG_H
#define RT_RNG_H

#include <stdint.h>

enum {
    RT_SITE_ATTEN    = 0x11u,
    RT_SITE_MATERIAL = 0x22u,
    RT_SITE_SAMPLING = 0x33u,
    RT_SITE_IMAGE    = 0x44u
};

/* splitmix64 finalizer: a bijective 64-bit mix. */
uint64_t rt_rng_mix64(uint64_t z);

/* Domain-separated key: site, work-item index, and per-item stream index
   (e.g. sample number) map to a unique draw. */
uint64_t rt_rng_key(unsigned site, uint64_t item, uint64_t stream);

/* 53-bit uniform in [0,1) and the raw 64-bit draw behind it. */
double   rt_rng_unit(uint64_t key);
uint64_t rt_rng_u64(uint64_t key);

#endif /* RT_RNG_H */
