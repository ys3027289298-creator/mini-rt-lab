#!/bin/sh
# negative_controls.sh -- prove the acceptance checks are not vacuous.
#
# Run from kernel/:   ./scripts/negative_controls.sh
#
# Each control injects one defect into a throwaway copy of the unified source
# and asserts that the corresponding check FAILS. A suite that passes both the
# real source and every mutant is not measuring anything; this script requires
# the opposite: real source PASS, mutants FAIL.
#
# The working tree is never modified; mutants are built in build/neg/.
set -eu
cd "$(dirname "$0")/.."
BUILD=build/neg
mkdir -p "$BUILD"

say() { printf '%s\n' "$*"; }
verdict=0
expect_fail() { # expect_fail <label> <command...>: pass if the command fails
    label=$1; shift
    if "$@" > "$BUILD/$label.log" 2>&1; then
        printf 'FAIL  control %s: defect was NOT detected\n' "$label"
        verdict=1
    else
        printf 'PASS  control %s: defect detected\n' "$label"
    fi
}

# --- control A: make the unified source stateful (a shared counter) ---------
# A counter-based source with a mutable draw counter is order-dependent; the
# schedule check must catch it.
cat > "$BUILD/rt_unified_stateful.c" <<'EOF'
#include <math.h>
#include "rt_kernel.h"
#include "rt_rng.h"
#include "rt_shim.h"
static uint64_t g_counter;              /* injected mutable state */
uint64_t rt_rng_configure(uint64_t s) { g_counter = s; return s; }
RtAttenSample rt_atten_sample(double sigma, uint64_t key) {
    (void)key;
    RtAttenSample r;
    double u = rt_rng_unit(g_counter++);
    r.t = -log1p(-u) / sigma;
    r.weight = exp(-sigma * r.t);
    return r;
}
double rt_material_mean_cos(int n, uint64_t key) {
    (void)key; if (n <= 0) return 0.0;
    double acc = 0.0; float v[3];
    for (int i = 0; i < n; i++) {
        float u1 = (float)rt_rng_unit(g_counter++);
        float u2 = (float)rt_rng_unit(g_counter++);
        rt_shim_cosine_sample(u1, u2, v); acc += (double)v[2];
    }
    return acc / (double)n;
}
double rt_strat_quad(int strata, uint64_t key) {
    (void)key; if (strata <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < strata; i++)
        sum += rt_shim_f_quad(((double)i + rt_rng_unit(g_counter++)) / (double)strata);
    return sum / (double)strata;
}
EOF
cc -O2 -std=c11 -ffp-contract=off -fno-fast-math -Iinclude \
   -DPROBE_SITE=0 -DPROBE_ITEMS=256 \
   tools/site_probe.c src/rt_rng.c "$BUILD/rt_unified_stateful.c" \
   src/rt_shim_material.c src/rt_shim_sample.c \
   -o "$BUILD/stateful_site0" -lm -pthread 2>/dev/null
# The unified suite's schedule check must reject this build.  Reproduce that
# check here directly.
a=$("$BUILD/stateful_site0" 1 256 seq | sha256sum | awk '{print $1}')
b=$("$BUILD/stateful_site0" 3 256 strided | sha256sum | awk '{print $1}')
if [ "$a" != "$b" ]; then
    printf 'PASS  control stateful_order_dependence: defect detected\n'
else
    printf 'FAIL  control stateful_order_dependence: defect was NOT detected\n'; verdict=1
fi

# --- control B: replace the mixer with a low-quality one -------------------
# x ^= x << 1 is affine over GF(2); its low bits are almost constant, so the
# uniformity / bit-balance tests must catch it.
cat > "$BUILD/rt_rng_weak.c" <<'EOF'
#include <stdint.h>
#include "rt_rng.h"
uint64_t rt_rng_mix64(uint64_t z) { return z ^ (z << 1); }
uint64_t rt_rng_key(unsigned site, uint64_t item, uint64_t stream) {
    return (uint64_t)site * 1000003ULL + item * 97ULL + stream;
}
uint64_t rt_rng_u64(uint64_t key) { return rt_rng_mix64(key); }
double rt_rng_unit(uint64_t key) {
    return (double)(rt_rng_u64(key) >> 11) * (1.0 / 9007199254740992.0);
}
EOF
cc -O2 -std=c11 -ffp-contract=off -fno-fast-math -Iinclude \
   tools/rng_selftest.c "$BUILD/rt_rng_weak.c" \
   -o "$BUILD/rng_weak" -lm 2>/dev/null
expect_fail weak_mixer "$BUILD/rng_weak"

# --- control C: change a public signature ---------------------------------
# Compile sigcheck against a mutated copy of the header (never the real one).
# The signature mismatch must be a compile error.
cp -r include "$BUILD/include_bad"
sed -i 's/double rt_strat_quad(int strata, uint64_t key);/double rt_strat_quad(unsigned strata, uint64_t key);/' \
    "$BUILD/include_bad/rt_kernel.h"
if cmp -s include/rt_kernel.h "$BUILD/include_bad/rt_kernel.h"; then
    printf 'FAIL  control bad_signature: mutation did not apply\n'; verdict=1
else
    if cc -O2 -std=c11 -ffp-contract=off -fno-fast-math \
          -Werror=incompatible-pointer-types -I"$BUILD/include_bad" \
          -c tools/sigcheck.c -o "$BUILD/sigcheck_bad.o" 2> "$BUILD/sigcheck_bad.log"; then
        printf 'FAIL  control bad_signature: defect was NOT detected\n'; verdict=1
    else
        printf 'PASS  control bad_signature: defect detected\n'
    fi
fi

# --- control D: make one site's estimator biased ---------------------------
# Drop the jitter's stratification (a classic off-by-one), so the sampling
# site's mean moves away from 1/3 by more than 5 standard errors.
cat > "$BUILD/rt_unified_biased.c" <<'EOF'
#include <math.h>
#include "rt_kernel.h"
#include "rt_rng.h"
#include "rt_shim.h"
uint64_t rt_rng_configure(uint64_t s) { (void)s; return 0; }
RtAttenSample rt_atten_sample(double sigma, uint64_t key) {
    RtAttenSample r; double u = rt_rng_unit(key);
    r.t = -log1p(-u) / sigma; r.weight = exp(-sigma * r.t); return r;
}
double rt_material_mean_cos(int n, uint64_t key) {
    if (n <= 0) return 0.0; double acc = 0.0; float v[3];
    for (int i = 0; i < n; i++) {
        float u1 = (float)rt_rng_unit(rt_rng_key(RT_SITE_MATERIAL, key, (uint64_t)i * 2u));
        float u2 = (float)rt_rng_unit(rt_rng_key(RT_SITE_MATERIAL, key, (uint64_t)i * 2u + 1u));
        rt_shim_cosine_sample(u1, u2, v); acc += (double)v[2];
    }
    return acc / (double)n;
}
double rt_strat_quad(int strata, uint64_t key) {
    if (strata <= 0) return 0.0; double sum = 0.0;
    for (int i = 0; i < strata; i++) {
        double u = ((double)i + 0.5 * rt_rng_unit(rt_rng_key(RT_SITE_SAMPLING, key, (uint64_t)i)))
                   / (double)strata;              /* injected bias: jitter halved */
        sum += rt_shim_f_quad(u);
    }
    return sum / (double)strata;
}
EOF
cc -O2 -std=c11 -ffp-contract=off -fno-fast-math -Iinclude \
   tools/site_stats.c src/rt_rng.c "$BUILD/rt_unified_biased.c" \
   src/rt_shim_material.c src/rt_shim_sample.c \
   -o "$BUILD/biased_stats" -lm -pthread 2>/dev/null
expect_fail biased_estimator "$BUILD/biased_stats" 2 16 20000

say ""
if [ "$verdict" = 0 ]; then
    say "RESULT: PASS (every injected defect was detected)"
    exit 0
else
    say "RESULT: FAIL (a defect slipped through)"
    exit 1
fi
