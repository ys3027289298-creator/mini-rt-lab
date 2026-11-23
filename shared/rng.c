/*
 * Unified PCG32 random source — see rng.h.
 * Fixed-width integer arithmetic only: byte-identical streams everywhere.
 */
#include "rng.h"

void rng_seed(Rng* r, uint64_t seed, uint64_t stream) {
    /* Distinct odd increment per stream, as prescribed by the PCG paper. */
    r->inc = (stream << 1u) | 1u;
    r->state = 0u;
    (void)rng_next_u32(r);
    r->state += seed;
    (void)rng_next_u32(r);
}

uint32_t rng_next_u32(Rng* r) {
    uint64_t oldstate = r->state;
    r->state = oldstate * 6364136223846793005ULL + r->inc;
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << (((uint32_t)(0) - rot) & 31u));
}

float rng_next_float(Rng* r) {
    /* 24 mantissa bits -> exactly representable, division by 2^24. */
    uint32_t bits = rng_next_u32(r) >> 8u;
    return (float)bits / 16777216.0f;
}

double rng_next_double(Rng* r) {
    /* 53 mantissa bits from two 32-bit draws, division by 2^53. */
    uint64_t hi = (uint64_t)(rng_next_u32(r) >> 5u); /* 27 bits */
    uint64_t lo = (uint64_t)(rng_next_u32(r) >> 6u); /* 26 bits */
    uint64_t bits = (hi << 26u) | lo;                /* 53 bits */
    return (double)bits / 9007199254740992.0;
}
