#ifndef PERF_RNG_H
#define PERF_RNG_H

/*
 * Deterministic PCG-XSH-RR 64/32 source for the performance harness.
 *
 * Fixed-width integer arithmetic only: no libc rand(), no clock.  The same
 * (seed, stream) pair yields the same sequence on every run and every
 * compiler, so a tier's image hash is evidence about the workload, not about
 * wall-clock luck.
 *
 * Header-only and static inline on purpose: each shim links its own copy and
 * there is no shared mutable state between the three rendering paths.
 */
#include <stdint.h>

typedef struct {
    uint64_t state;
    uint64_t inc;
} PerfRng;

static inline uint32_t perf_rng_next(PerfRng* r) {
    uint64_t old = r->state;
    r->state = old * UINT64_C(6364136223846793005) + r->inc;
    uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = (uint32_t)(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((0u - rot) & 31u));
}

static inline void perf_rng_seed(PerfRng* r, uint64_t seed, uint64_t stream) {
    r->inc = (stream << 1u) | 1u;
    r->state = 0u;
    (void)perf_rng_next(r);
    r->state += seed;
    (void)perf_rng_next(r);
}

/* Uniform float in [0, 1), 24 bits of entropy (never rounds up to 1). */
static inline float perf_rng_float(PerfRng* r) {
    return (float)(perf_rng_next(r) >> 8u) * (1.0f / 16777216.0f);
}

#endif /* PERF_RNG_H */
