/* rng_selftest.c -- properties of the unified random source. */
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include "rt_rng.h"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* Purity: the same key gives the same draw, every time, no state. */
    int pure = 1;
    for (uint64_t k = 0; k < 4096; k++) {
        uint64_t a = rt_rng_u64(rt_rng_key(RT_SITE_ATTEN, k, 0));
        uint64_t b = rt_rng_u64(rt_rng_key(RT_SITE_ATTEN, k, 0));
        if (a != b) pure = 0;
    }
    check("purity_same_key_same_draw", pure);

    /* Domain separation: keys from distinct (site, item, stream) never collide
       over a dense sweep. */
    int distinct = 1;
    uint64_t seen[4096];
    int n = 0;
    for (unsigned site = 0; site < 4; site++)
        for (uint64_t item = 0; item < 32; item++)
            for (uint64_t stream = 0; stream < 32; stream++) {
                uint64_t k = rt_rng_key(site, item, stream);
                for (int j = 0; j < n; j++) if (seen[j] == k) distinct = 0;
                seen[n++] = k;
            }
    check("domain_separation_distinct_keys", distinct);

    /* Uniformity: 53-bit values fill [0,1), never 1.0, mean ~ 0.5. */
    int in_range = 1;
    double sum = 0.0;
    const int m = 200000;
    for (int i = 0; i < m; i++) {
        double u = rt_rng_unit(rt_rng_key(RT_SITE_IMAGE, (uint64_t)i, 0));
        if (!(u >= 0.0 && u < 1.0)) in_range = 0;
        sum += u;
    }
    check("uniform_in_half_open_unit", in_range);
    check("uniform_mean", fabs(sum / m - 0.5) < 0.01);

    /* Bit balance: every bit position is set in roughly half the draws. */
    int bits_ok = 1;
    for (int bit = 0; bit < 64; bit++) {
        int ones = 0;
        for (int i = 0; i < 4096; i++)
            ones += (int)((rt_rng_u64(rt_rng_key(RT_SITE_IMAGE, (uint64_t)i, 1)) >> bit) & 1u);
        if (ones < 4096 / 3 || ones > 4096 * 2 / 3) bits_ok = 0;
    }
    check("bit_balance", bits_ok);

    /* No short cycle: a 4096-step walk revisits no value. */
    int no_repeat = 1;
    for (int i = 0; i < 4096; i++) {
        uint64_t a = rt_rng_u64((uint64_t)i);
        for (int j = 0; j < i; j++) if (rt_rng_u64((uint64_t)j) == a) no_repeat = 0;
    }
    check("no_short_cycle", no_repeat);

    /* Layout independence: interleaving draws from other keys cannot change a
       given key's value (there is no state to advance). */
    uint64_t before = rt_rng_u64(rt_rng_key(RT_SITE_SAMPLING, 7, 3));
    for (int i = 0; i < 1000; i++) (void)rt_rng_u64(rt_rng_key(RT_SITE_MATERIAL, (uint64_t)i, 0));
    uint64_t after = rt_rng_u64(rt_rng_key(RT_SITE_SAMPLING, 7, 3));
    check("draw_order_independent", before == after);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
