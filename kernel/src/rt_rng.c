/* rt_rng.c -- the single random source.  See include/rt_rng.h. */
#include "rt_rng.h"

uint64_t rt_rng_mix64(uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

uint64_t rt_rng_key(unsigned site, uint64_t item, uint64_t stream) {
    /* Fold each field through the mixer so adjacent items, streams and sites
       cannot alias, then combine. */
    uint64_t a = rt_rng_mix64((uint64_t)site * 0x100000001b3ULL + 0x9e3779b9u);
    uint64_t b = rt_rng_mix64(item);
    uint64_t c = rt_rng_mix64(stream ^ 0xd1b54a32d192ed03ULL);
    return rt_rng_mix64(a ^ (b + 0x165667b19e3779f9ULL) ^ (c << 1));
}

uint64_t rt_rng_u64(uint64_t key) { return rt_rng_mix64(key); }

double rt_rng_unit(uint64_t key) {
    /* 53 significant bits: the largest value is 1 - 2^-53, so it cannot
       round to 1.0 in double. */
    return (double)(rt_rng_u64(key) >> 11) * (1.0 / 9007199254740992.0);
}
