
/* sample: stratified and low-discrepancy sampling with variance estimate. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef RT_UNIFIED_RNG
#include "rng.h"
#include "scene.h"

/* Unified stream for the sampling module. Off-switch compiles this out and
   keeps every original code path untouched. */
static Rng g_rng;
#endif

#define N_MAX 256

typedef struct { double mean; double var; int n; } Stat;

static unsigned radical_inverse(unsigned i, unsigned base) {
    unsigned rev = 0, base_n = base;
    while (i > 0) { rev += (i % base) * (base_n / base); i /= base; base_n *= base; }
    return rev;
}

static double halton(unsigned i, unsigned base) { return (double)radical_inverse(i, base) / 4294967296.0 * base; }
static double halton2(unsigned i, unsigned base) {
    double f = 1.0, r = 0.0; unsigned k = i;
    while (k > 0) { f /= base; r += f * (k % base); k /= base; }
    return r;
}

/* estimate integral of f over [0,1) with stratified samples */
static Stat stratify(int strata, double (*f)(double)) {
    Stat s; memset(&s, 0, sizeof(s));
    double sum = 0.0, sum2 = 0.0;
    int n = 0;
    for (int i = 0; i < strata; i++) {
        double u = ((double)i + 0.5) / (double)strata;
#ifdef RT_UNIFIED_RNG
        /* jitter inside the stratum from the unified stream */
        u = ((double)i + rng_next_double(&g_rng)) / (double)strata;
#endif
        double v = f(u);
        sum += v; sum2 += v * v; n++;
    }
    s.n = n;
    s.mean = sum / n;
    double var = (n > 1) ? (sum2 / n - s.mean * s.mean) * n / (n - 1) : 0.0;
    s.var = var > 0 ? var : 0;
    return s;
}

static double f_quad(double u) { return u * u; }
static double f_const(double u) { (void)u; return 1.0; }

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
#ifdef RT_UNIFIED_RNG
    rng_seed(&g_rng, SCENE_SEED, STREAM_SAMPLE);
#endif
    check("halton_in_unit", halton2(1, 2) >= 0.0 && halton2(1, 2) < 1.0);
    check("halton_first", fabs(halton2(0, 2)) < 1e-12);

    Stat a = stratify(64, f_quad);
    check("quad_mean_close", fabs(a.mean - 1.0 / 3.0) < 1e-3);

    Stat b = stratify(16, f_quad);
    Stat c = stratify(64, f_quad);
    check("variance_drops", c.var <= b.var + 1e-6);

    Stat d = stratify(32, f_const);
    check("const_exact", fabs(d.mean - 1.0) < 1e-12 && d.var < 1e-12);

    check("halton_distinct", halton2(1, 2) != halton2(2, 2));
    check("halton_monotone_base", halton2(1, 2) > halton2(0, 2));

#ifdef RT_UNIFIED_RNG
    /* Unified stream must be reproducible: same (seed, stream) -> same draws.
       Checked on two fresh RNGs so g_rng stays untouched for stratify(). */
    Rng ra, rb;
    rng_seed(&ra, SCENE_SEED, STREAM_SAMPLE);
    rng_seed(&rb, SCENE_SEED, STREAM_SAMPLE);
    int rep_ok = 1;
    for (int i = 0; i < 4096; i++) {
        if (rng_next_u32(&ra) != rng_next_u32(&rb)) { rep_ok = 0; break; }
    }
    check("rng_reproducible", rep_ok);

    Rng rs, rm;
    rng_seed(&rs, SCENE_SEED, STREAM_SAMPLE);
    rng_seed(&rm, SCENE_SEED, STREAM_MATERIAL);
    check("rng_streams_distinct", rng_next_u32(&rs) != rng_next_u32(&rm));
    printf("rng=unified\n");
#endif

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
