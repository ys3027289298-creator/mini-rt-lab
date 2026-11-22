#!/usr/bin/env bash
#
# consistency_test.sh — temporal consistency gate for the renderer.
#
# Requirements under test:
#   R1  Same parameters, repeated renders -> zero per-pixel/per-channel diff
#   R2  Rendering with different thread counts -> byte-identical output,
#       including thread counts that do not divide the frame height
#   R3  When frames differ, compare reports the FIRST differing pixel with
#       its (x, y) and channel (via --selftest and a fault-injected run)
#   R4  Fault-injection mode (--rng=global, one shared clock-seeded RNG)
#       must actually be detected as non-deterministic, proving the harness
#       can catch the "global state" bug class.
#
# Exit status is non-zero if any check fails.
set -u

cd "$(dirname "$0")/.."

W=64
H=48
SPP=16
SEED=7
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

failures=0
checks=0

pass() { checks=$((checks + 1)); printf 'PASS %s\n' "$1"; }
fail() { checks=$((checks + 1)); failures=$((failures + 1)); printf 'FAIL %s\n' "$1"; }

# expect_zero LABEL -- command...
expect_zero() {
    local label="$1"; shift
    if "$@"; then pass "$label"; else fail "$label"; fi
}

# expect_nonzero LABEL -- command...
expect_nonzero() {
    local label="$1"; shift
    if "$@" >/dev/null 2>&1; then fail "$label"; else pass "$label"; fi
}

echo "== R1: same parameters, two independent runs =="
./render -w "$W" -h "$H" -s "$SPP" -t 4 --seed "$SEED" -o "$TMP/a.pfm" 2>/dev/null
./render -w "$W" -h "$H" -s "$SPP" -t 4 --seed "$SEED" -o "$TMP/b.pfm" 2>/dev/null
expect_zero "repeat_render_identical" ./compare "$TMP/a.pfm" "$TMP/b.pfm"
# Independent byte-level cross-check, not relying on compare itself.
if cmp -s "$TMP/a.pfm" "$TMP/b.pfm"; then
    pass "repeat_render_bytes_equal"
else
    fail "repeat_render_bytes_equal"
fi

echo "== R2: thread count changes must not change any pixel =="
# 48 rows: t=1,2,4,8 divide evenly; 3 and 7 do not; 64 > h creates idle bands.
for t in 1 2 3 4 7 8 64; do
    ./render -w "$W" -h "$H" -s "$SPP" -t "$t" --seed "$SEED" \
        -o "$TMP/t$t.pfm" 2>/dev/null
    expect_zero "threads_${t}_matches_baseline" \
        ./compare "$TMP/a.pfm" "$TMP/t$t.pfm"
done
# And confirm the reference PFM itself is stable under repeated rebuilds.
sha1=$(sha256sum "$TMP/a.pfm" | cut -d' ' -f1)
sha2=$(sha256sum "$TMP/t7.pfm" | cut -d' ' -f1)
if [ "$sha1" = "$sha2" ]; then
    pass "sha256_identical_across_thread_counts"
else
    fail "sha256_identical_across_thread_counts"
fi

echo "== R3: mismatch must report the first differing pixel and channel =="
expect_zero "compare_selftest" ./compare --selftest
# Dimension mismatch is a usage-level error (exit 2), never a silent "equal".
./render -w 32 -h 24 -s "$SPP" -t 2 --seed "$SEED" -o "$TMP/small.pfm" 2>/dev/null
expect_nonzero "dimension_mismatch_detected" \
    ./compare "$TMP/a.pfm" "$TMP/small.pfm"

echo "== R4: fault injection (global clock-seeded RNG) must be caught =="
./render -w "$W" -h "$H" -s "$SPP" -t 4 --rng=global -o "$TMP/g1.pfm" 2>/dev/null
./render -w "$W" -h "$H" -s "$SPP" -t 4 --rng=global -o "$TMP/g2.pfm" 2>/dev/null
if ./compare "$TMP/g1.pfm" "$TMP/g2.pfm"; then
    fail "global_rng_nondeterminism_detected"
    echo "  (global RNG runs compared equal — fault injection not working)"
else
    rc=$?
    if [ "$rc" -eq 1 ]; then
        pass "global_rng_nondeterminism_detected"
    else
        fail "global_rng_nondeterminism_detected (unexpected exit $rc)"
    fi
fi
# Scheduling-dependent path: same global run vs. the deterministic frame.
expect_nonzero "global_vs_deterministic_differs" \
    ./compare "$TMP/a.pfm" "$TMP/g1.pfm"
# And with a different thread count the global path must diverge too.
./render -w "$W" -h "$H" -s "$SPP" -t 7 --rng=global -o "$TMP/g3.pfm" 2>/dev/null
expect_nonzero "global_rng_thread_sensitive" \
    ./compare "$TMP/g1.pfm" "$TMP/g3.pfm"

echo
echo "checks=$checks failures=$failures"
if [ "$failures" -ne 0 ]; then
    echo "CONSISTENCY TEST FAILED"
    exit 1
fi
echo "CONSISTENCY TEST PASSED"
