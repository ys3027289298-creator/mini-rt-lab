#ifndef MINI_RT_RNG_H
#define MINI_RT_RNG_H

/*
 * Unified random source for sampling, material and geometry.
 *
 * PCG-XSH-RR 64/32 (Melissa O'Neill's PCG): every operation uses fixed-width
 * integer types, so the produced sequence is identical on every machine and
 * every compiler that implements stdint correctly. No libc rand(), no clock,
 * no platform-dependent state.
 *
 * Each consumer seeds an independent stream from the same scene seed; the
 * stream id makes the three sequences reproducible and disjoint.
 */

#include <stdint.h>

typedef struct {
    uint64_t state;
    uint64_t inc;
} Rng;

/* Seed a stream. Two Rng seeded with the same (seed, stream) pair emit the
 * exact same sequence forever. */
void rng_seed(Rng* r, uint64_t seed, uint64_t stream);

/* Raw 32-bit draw. */
uint32_t rng_next_u32(Rng* r);

/* Uniform float in [0, 1), 24 bits of entropy. */
float rng_next_float(Rng* r);

/* Uniform double in [0, 1), 53 bits of entropy. */
double rng_next_double(Rng* r);

#endif /* MINI_RT_RNG_H */
