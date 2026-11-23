/* sigcheck.c -- compile-time proof that the public signatures are unchanged.
 *
 * Each assignment only compiles if the declaration has exactly the expected
 * type, so a changed signature, parameter type or return type is a build
 * error.  There is no runtime logic here; returning 0 means it compiled.
 */
#include <stdio.h>
#include "rt_kernel.h"

typedef RtAttenSample (*atten_fn)(double, uint64_t);
typedef double         (*mean_cos_fn)(int, uint64_t);
typedef double         (*strat_fn)(int, uint64_t);
typedef uint64_t       (*configure_fn)(uint64_t);

static atten_fn      const p_atten     = rt_atten_sample;
static mean_cos_fn   const p_mean_cos  = rt_material_mean_cos;
static strat_fn      const p_strat     = rt_strat_quad;
static configure_fn  const p_configure = rt_rng_configure;

static RtAttenSample* const check_record = (RtAttenSample*)0;

int main(void) {
    printf("sig_atten=%d\n", p_atten != NULL);
    printf("sig_mean_cos=%d\n", p_mean_cos != NULL);
    printf("sig_strat=%d\n", p_strat != NULL);
    printf("sig_configure=%d\n", p_configure != NULL);
    printf("record_fields_t=%d\n", (int)sizeof(RtAttenSample) == (int)(2 * sizeof(double)));
    printf("record_ptr_t=%d\n", check_record == NULL);
    printf("all_ok=1\n");
    return 0;
}
