#!/bin/sh
# portability.sh -- results must not depend on compiler, optimisation or the
# machine that built them.
#
# Run from kernel/:   ./scripts/portability.sh
#
# Builds the unified scene fingerprint with gcc and clang at -O0/-O1/-O2/-O3
# and requires every binary to print the identical scene_hash.  A generator
# whose values depend on draw order, word size quirks or floating-point
# contraction fails this; the counter-based source passes it because each
# value is a pure function of its key.
#
# It also checks the OFF switch across compilers: the legacy build's
# single-thread output must equal the frozen baseline regardless of who
# compiled it.
set -eu
cd "$(dirname "$0")/.."
BUILD=build/portability
EVID=evidence
mkdir -p "$BUILD" "$EVID"

compilers=""
command -v gcc   >/dev/null 2>&1 && compilers="$compilers gcc"
command -v clang >/dev/null 2>&1 && compilers="$compilers clang"
[ -n "$compilers" ] || { echo "FAIL no compiler found"; exit 1; }

fail=0
hashes=""
: > "$EVID/portability.txt"

for cc in $compilers; do
    for opt in -O0 -O1 -O2 -O3; do
        out="$BUILD/scene_${cc}${opt}"
        $cc $opt -std=c11 -ffp-contract=off -fno-fast-math -Iinclude \
            tools/scene_hash.c src/rt_rng.c src/rt_unified.c \
            src/rt_shim_material.c src/rt_shim_sample.c \
            -o "$out" -lm -pthread 2>/dev/null
        line=$("$out" 1 256 seq)
        printf '%s\n' "$line" >> "$EVID/portability.txt"
        hashes="$hashes $(printf '%s' "$line" | awk '{print $2}')"
    done
done

uniq_hashes=$(printf '%s\n' $hashes | sort -u | wc -l)
nbuilds=$(printf '%s\n' $hashes | wc -w)
if [ "$uniq_hashes" = 1 ]; then
    printf 'PASS  unified scene_hash identical across %s compiler/opt builds\n' "$nbuilds"
else
    printf 'FAIL  unified scene_hash differs across builds (see %s)\n' "$EVID/portability.txt"
    fail=1
fi

# OFF switch across compilers
for cc in $compilers; do
    for site in 0 1 2; do
        # rebuild per site
        $cc -O2 -std=c11 -ffp-contract=off -fno-fast-math -Iinclude \
            -DPROBE_SITE=$site tools/site_probe.c src/rt_rng.c src/rt_legacy.c \
            src/rt_shim_material.c src/rt_shim_sample.c \
            -o "$BUILD/legacy_${cc}_s$site" -lm -pthread 2>/dev/null
        if "$BUILD/legacy_${cc}_s$site" 1 256 seq | cmp -s - \
           "baseline/legacy/site$site.threads1.txt"; then
            printf 'PASS  legacy site%s baseline reproduced with %s\n' "$site" "$cc"
        else
            printf 'FAIL  legacy site%s baseline differs with %s\n' "$site" "$cc"
            fail=1
        fi
    done
done

[ "$fail" = 0 ] && { echo "RESULT: PASS"; exit 0; } || { echo "RESULT: FAIL"; exit 1; }
