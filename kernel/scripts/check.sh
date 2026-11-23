#!/bin/sh
# check.sh -- one-command acceptance run for the unified random source.
#
# Run from kernel/:   ./scripts/check.sh        (or: make test)
#
# Verifies, in order:
#   1. the three existing modules still build and their self-test stdout is
#      byte-identical to the frozen baseline;
#   2. the public API is unchanged (header hash, compile-time signature
#      assertions, exported symbol table identical between both builds);
#   3. the unified source passes its own property tests;
#   4. OFF switch: the legacy ("unified source disabled") build reproduces the
#      frozen pre-change goldens byte-for-byte, for all three sites;
#   5. unified path: every site is thread-count independent (1..8 threads,
#      including non-divisors) and repeat-stable;
#   6. error data for all three sites, before vs after, against analytic
#      targets, each within a 5-standard-error bound;
#   7. no mutable generator state leaks into the unified objects;
#   8. ASan/UBSan on both builds, and TSan on the unified threaded probe.
#
# Exit status 0 only if every check passes.  Nothing tracked is rewritten:
# generated files land in build/ and evidence/.
set -eu

cd "$(dirname "$0")/.."
ROOT=$(pwd)
EVID="$ROOT/evidence"
BUILD="$ROOT/build"
SRCROOT="$ROOT/.."
mkdir -p "$EVID" "$BUILD"

fail=0
say()  { printf '%s\n' "$*"; }
ok()   { printf 'PASS  %s\n' "$*"; }
bad()  { printf 'FAIL  %s\n' "$*"; fail=1; }
same() { if cmp -s "$2" "$3"; then ok "$1"; else bad "$1"; fi; }

say "=== kernel acceptance $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="

# ---------------------------------------------------------------- 1. modules
say "--- 1. existing modules unchanged ---"
for m in bvh sample material; do
    make -s -C "$SRCROOT/$m" clean >/dev/null 2>&1 || true
    make -s -C "$SRCROOT/$m" >/dev/null 2>&1
    "$SRCROOT/$m/$m" > "$EVID/module_$m.txt" 2>/dev/null || true
    same "module_$m selftest byte-identical" \
         "$ROOT/baseline/module_tests/$m.txt" "$EVID/module_$m.txt"
done
git -C "$SRCROOT" diff --stat -- bvh sample material > "$EVID/module_git_diff.txt" 2>/dev/null || true
if [ -s "$EVID/module_git_diff.txt" ]; then
    bad "bvh/ sample/ material/ tracked files unmodified"
else
    ok "bvh/ sample/ material/ tracked files unmodified"
fi

# ------------------------------------------------------------- 2. public API
say "--- 2. public API frozen ---"
if [ "$(sha256sum "$ROOT/include/rt_kernel.h" | awk '{print $1}')" = \
     "$(awk '{print $1}' "$ROOT/baseline/public_api.sha256")" ]; then
    ok "rt_kernel.h sha256 unchanged"
else
    bad "rt_kernel.h sha256 unchanged"
fi

if "$BUILD/sigcheck" > "$EVID/sigcheck.txt"; then
    ok "compile-time signature assertions"
else
    bad "compile-time signature assertions"
fi

if diff -q "$BUILD/symbols_legacy.txt" "$BUILD/symbols_unified.txt" >/dev/null; then
    ok "exported symbol tables identical (legacy vs unified)"
else
    bad "exported symbol tables identical (legacy vs unified)"
fi
cp "$BUILD/symbols_unified.txt" "$EVID/symbols.txt"

# ------------------------------------------------------------ 3. RNG itself
say "--- 3. unified source property tests ---"
if "$BUILD/rng_selftest" > "$EVID/rng_selftest.txt"; then
    ok "rng_selftest"
else
    bad "rng_selftest"
fi

# ------------------------------------------------------- 4. off switch (gold)
say "--- 4. OFF switch reproduces frozen pre-change output ---"
for site in 0 1 2; do
    "$BUILD/legacy_site$site" 1 256 > "$EVID/legacy_site$site.txt"
    same "site$site legacy == frozen baseline" \
         "$ROOT/baseline/legacy/site$site.threads1.txt" "$EVID/legacy_site$site.txt"
done

# ---------------------------------------------------- 5. thread independence
say "--- 5. unified thread-count independence ---"
for site in 0 1 2; do
    "$BUILD/unified_site$site" 1 256 > "$BUILD/u${site}_t1.txt"
    site_ok=1
    for t in 2 3 4 5 7 8; do
        "$BUILD/unified_site$site" "$t" 256 > "$BUILD/u${site}_t$t.txt"
        cmp -s "$BUILD/u${site}_t1.txt" "$BUILD/u${site}_t$t.txt" || site_ok=0
    done
    "$BUILD/unified_site$site" 7 256 > "$BUILD/u${site}_t7b.txt"
    cmp -s "$BUILD/u${site}_t7.txt" "$BUILD/u${site}_t7b.txt" || site_ok=0
    if [ "$site_ok" = 1 ]; then
        ok "site$site identical at 1,2,3,4,5,7,8 threads and on repeat"
    else
        bad "site$site identical across thread counts"
    fi
    cp "$BUILD/u${site}_t1.txt" "$EVID/unified_site$site.txt"
done

say "--- 5b. legacy path is schedule-dependent (defect contrast) ---"
# Race-free emulation: same work items, different visit order.  No threads are
# involved, so the difference can only come from generator state.
for site in 0 1 2; do
    h1=$("$BUILD/legacy_site$site" 1 256 seq | sha256sum | awk '{print $1}')
    h2=$("$BUILD/legacy_site$site" 3 256 strided | sha256sum | awk '{print $1}')
    if [ "$h1" = "$h2" ]; then
        bad "site$site legacy varies with visit order (defect not reproduced)"
    else
        ok "site$site legacy varies with visit order (defect reproduced)"
    fi
    u1=$("$BUILD/unified_site$site" 1 256 seq | sha256sum | awk '{print $1}')
    u2=$("$BUILD/unified_site$site" 3 256 strided | sha256sum | awk '{print $1}')
    if [ "$u1" = "$u2" ]; then
        ok "site$site unified invariant under visit order"
    else
        bad "site$site unified invariant under visit order"
    fi
done

# ---------------------------------------------- 5c. whole-scene fingerprint
say "--- 5c. scene fingerprint (all three sites at once) ---"
: > "$EVID/scene_hash.txt"
ref=$("$BUILD/scene_hash_unified" 1 256 seq | awk '{print $2}')
for cfg in "1 256 seq" "2 256 seq" "3 256 seq" "4 256 seq" "7 256 seq" \
           "8 256 seq" "2 256 strided" "3 256 strided" "8 256 strided"; do
    line=$("$BUILD/scene_hash_unified" $cfg)
    printf '%s\n' "$line" >> "$EVID/scene_hash.txt"
    h=$(printf '%s' "$line" | awk '{print $2}')
    if [ "$h" = "$ref" ]; then
        ok "scene_hash $cfg == $ref"
    else
        bad "scene_hash $cfg = $h, expected $ref"
    fi
done
{
    echo "-- legacy path, for contrast"
    "$BUILD/scene_hash_legacy" 1 256 seq
    "$BUILD/scene_hash_legacy" 3 256 strided
} >> "$EVID/scene_hash.txt"

# --------------------------------------------------------------- 6. accuracy
say "--- 6. error data, before vs after ---"
: > "$EVID/accuracy.txt"
acc_ok=1
run_pair() { # run_pair <site> <param> <label>
    site=$1; param=$2; label=$3
    "$BUILD/stats_legacy"  "$site" "$param" 20000 > "$BUILD/acc_${label}_legacy.txt"  || acc_ok=0
    "$BUILD/stats_unified" "$site" "$param" 20000 > "$BUILD/acc_${label}_unified.txt" || acc_ok=0
    {
        echo "### $label (site $site, param $param, 20000 items)"
        echo "-- before: legacy per-site RNG"
        cat "$BUILD/acc_${label}_legacy.txt"
        echo "-- after: unified counter-based source"
        cat "$BUILD/acc_${label}_unified.txt"
        echo
    } >> "$EVID/accuracy.txt"
}
run_pair 0 1.7 attenuation
run_pair 1 8   material
run_pair 2 16  sampling

# Convergence: at 4x the items the standard error must shrink, and bias must
# stay inside the same bound.  Printed, not gated.
{
    echo "### convergence (4x items)"
    for spec in "0 1.7" "1 8" "2 16"; do
        set -- $spec
        printf '80k items: '
        "$BUILD/stats_unified" "$1" "$2" 80000 | tr '\n' ' '
        echo
    done
    echo
} >> "$EVID/accuracy.txt"
if [ "$acc_ok" = 1 ]; then
    ok "all three sites within 5 stderr of analytic target (both paths)"
else
    bad "accuracy bound"
fi

# ----------------------------------------------- 7. no mutable global state
say "--- 7. unified objects carry no mutable generator state ---"
if nm "$BUILD/rt_impl_u.o" | grep -qE ' [bBdD] '; then
    bad "rt_unified.o has writable data symbols"
else
    ok "rt_unified.o has no writable data symbols"
fi
if nm "$BUILD/rt_rng_u.o" | grep -qE ' [bBdD] '; then
    bad "rt_rng.o has writable data symbols"
else
    ok "rt_rng.o has no writable data symbols"
fi
# record the contrast (legacy does have state)
{
    echo "--- writable data symbols per object ---"
    for o in rt_impl_l.o rt_impl_u.o rt_rng_l.o rt_rng_u.o; do
        printf '%-14s ' "$o"
        nm "$BUILD/$o" 2>/dev/null | awk '$2=="b"||$2=="B"||$2=="d"||$2=="D"{printf "%s(%s) ", $3, $2} END{print ""}'
    done
} > "$EVID/global_state.txt"

# ------------------------------------------------------------ 8. sanitizers
say "--- 8. sanitizers ---"
for cfg in legacy unified; do
    impl=$([ "$cfg" = legacy ] && echo src/rt_legacy.c || echo src/rt_unified.c)
    cc -O1 -g -std=c11 -ffp-contract=off -fno-fast-math \
       -fsanitize=address,undefined -fno-sanitize-recover=all \
       -Iinclude -DPROBE_SITE=1 -DPROBE_ITEMS=64 \
       tools/site_probe.c src/rt_rng.c "$impl" \
       src/rt_shim_material.c src/rt_shim_sample.c \
       -o "$BUILD/asan_$cfg" -lm -pthread 2>/dev/null
    if "$BUILD/asan_$cfg" 4 64 > /dev/null 2> "$BUILD/asan_${cfg}.log"; then
        ok "ASan/UBSan clean ($cfg)"
    else
        bad "ASan/UBSan ($cfg)"; sed -n '1,20p' "$BUILD/asan_${cfg}.log" || true
    fi
done

if command -v clang >/dev/null 2>&1; then
    if clang -O1 -g -std=c11 -ffp-contract=off -fno-fast-math -fsanitize=thread \
             -Iinclude -DPROBE_SITE=1 -DPROBE_ITEMS=256 \
             tools/site_probe.c src/rt_rng.c src/rt_unified.c \
             src/rt_shim_material.c src/rt_shim_sample.c \
             -o "$BUILD/tsan_unified" -lm -pthread 2>/dev/null && \
       "$BUILD/tsan_unified" 8 256 > /dev/null 2> "$BUILD/tsan.log"; then
        ok "TSan clean (unified, 8 threads)"
    else
        bad "TSan (unified, 8 threads)"; sed -n '1,20p' "$BUILD/tsan.log" 2>/dev/null || true
    fi
else
    say "note: clang not found, TSan check skipped"
fi

# ------------------------------------------------ 9. compiler portability
say "--- 9. compiler / optimisation portability ---"
if ./scripts/portability.sh > "$BUILD/portability.out" 2>&1; then
    ok "portability ($(grep -c '^PASS' "$BUILD/portability.out") checks)"
else
    bad "portability"; sed -n '1,20p' "$BUILD/portability.out"
fi

# ------------------------------------------ 10. controls are not vacuous
say "--- 10. negative controls (checks reject injected defects) ---"
if ./scripts/negative_controls.sh > "$BUILD/neg.out" 2>&1; then
    ok "all injected defects detected ($(grep -c '^PASS' "$BUILD/neg.out") controls)"
else
    bad "negative controls"; sed -n '1,20p' "$BUILD/neg.out"
fi

# ------------------------------------------------------------------ verdict
say ""
if [ "$fail" = 0 ]; then
    say "RESULT: PASS (all acceptance checks)"
    exit 0
else
    say "RESULT: FAIL (see FAIL lines above)"
    exit 1
fi
