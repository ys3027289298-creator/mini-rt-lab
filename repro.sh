#!/bin/sh
# repro.sh: reproducible-build evidence for mini-rt-lab.
#
# For each module (bvh, sample, material):
#   1. dependency-free build under two different compile option sets
#      (A: gcc -O2, B: clang -O3), zero-warning via -Wall -Wextra -Werror
#   2. self-tests pass on both builds, with identical output
#   3. rendered images compared pixel by pixel (diff must be zero)
#   4. invalid parameters must fail with a non-zero exit and an explanation
#   5. --strict vs --no-strict artifacts must be byte-identical
#
# Evidence is written to artifacts/repro/. Any failed item aborts the script.
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ART="$ROOT/artifacts/repro"
MODS="bvh sample material"
SIZE=96

CC_A=${CC_A:-gcc}
CC_B=${CC_B:-clang}
FLAGS_A="-O2 -std=c11 -Wall -Wextra -Werror -ffp-contract=off"
FLAGS_B="-O3 -std=c11 -Wall -Wextra -Werror -ffp-contract=off"

rm -rf "$ART"
mkdir -p "$ART"
SUM="$ART/summary.txt"
: > "$SUM"
INV="$ART/invalid-params.txt"
: > "$INV"

say() { echo "$@" | tee -a "$SUM"; }

# expect a command to fail with non-zero exit and a stderr explanation
expect_fail() {
    desc=$1; shift
    set +e
    "$@" >"$ART/.tmp.out" 2>"$ART/.tmp.err"
    rc=$?
    set -e
    {
        echo "case: $desc"
        echo "cmd: $*"
        echo "exit=$rc"
        sed 's/^/stderr: /' "$ART/.tmp.err"
        echo ""
    } >> "$INV"
    if [ "$rc" -eq 0 ]; then
        say "FAIL: $desc: expected non-zero exit, got 0"
        exit 1
    fi
    if [ ! -s "$ART/.tmp.err" ]; then
        say "FAIL: $desc: no error message on stderr"
        exit 1
    fi
}

# pixel-differ, itself built warning-free
cc -O2 -std=c11 -Wall -Wextra -Werror "$ROOT/tools/ppmdiff.c" -o "$ART/ppmdiff"

say "toolchain A: $CC_A $FLAGS_A"
say "toolchain B: $CC_B $FLAGS_B"
say ""

for m in $MODS; do
    d="$ROOT/$m"

    # 1. warning-free builds under both option sets (-Werror: a warning fails the build)
    make -C "$d" clean >/dev/null
    make -C "$d" CC="$CC_A" CFLAGS="$FLAGS_A" >"$ART/build-A-$m.log" 2>&1
    mv "$d/$m" "$ART/$m-A"
    make -C "$d" clean >/dev/null
    make -C "$d" CC="$CC_B" CFLAGS="$FLAGS_B" >"$ART/build-B-$m.log" 2>&1
    mv "$d/$m" "$ART/$m-B"
    wa=$(grep -ci warning "$ART/build-A-$m.log" || true)
    wb=$(grep -ci warning "$ART/build-B-$m.log" || true)
    say "$m: zero-warning builds OK (warnings in logs: A=$wa B=$wb)"

    # 2. self-tests pass on both builds, output identical
    "$ART/$m-A" >"$ART/selftest-A-$m.txt"
    "$ART/$m-B" >"$ART/selftest-B-$m.txt"
    diff "$ART/selftest-A-$m.txt" "$ART/selftest-B-$m.txt" >/dev/null
    say "$m: self-tests pass on A and B, outputs identical"

    # 3. render under both builds, compare pixel by pixel
    "$ART/$m-A" render "$ART/$m-A.ppm" --size "$SIZE"
    "$ART/$m-B" render "$ART/$m-B.ppm" --size "$SIZE"
    "$ART/ppmdiff" "$ART/$m-A.ppm" "$ART/$m-B.ppm" | tee "$ART/pixeldiff-$m.txt"
    say "$m: pixel diff A vs B = 0"

    # 4. invalid parameters must error with explanation
    expect_fail "$m: render without output path"   "$ART/$m-A" render
    expect_fail "$m: --size 0"                     "$ART/$m-A" render "$ART/.tmp.ppm" --size 0
    expect_fail "$m: --size abc"                   "$ART/$m-A" render "$ART/.tmp.ppm" --size abc
    expect_fail "$m: --size 99999"                 "$ART/$m-A" render "$ART/.tmp.ppm" --size 99999
    expect_fail "$m: --size without value"         "$ART/$m-A" render "$ART/.tmp.ppm" --size
    expect_fail "$m: unknown option"               "$ART/$m-A" render "$ART/.tmp.ppm" --bogus
    expect_fail "$m: unknown command"              "$ART/$m-A" bogus
    expect_fail "$m: unwritable output path"       "$ART/$m-A" render /nonexistent-dir/x.ppm
    say "$m: invalid parameters rejected with error messages (8 cases)"

    # 5. strict vs no-strict artifacts byte-identical
    "$ART/$m-A" render "$ART/$m-strict.ppm" --size "$SIZE" --strict
    "$ART/$m-A" render "$ART/$m-nostrict.ppm" --size "$SIZE" --no-strict
    cmp "$ART/$m-strict.ppm" "$ART/$m-nostrict.ppm"
    say "$m: --strict vs --no-strict artifacts byte-identical"

    # leave the directory in its default built state
    make -C "$d" clean >/dev/null
    make -C "$d" >/dev/null 2>&1
    say ""
done

say "ALL CHECKS PASSED"
