#!/bin/sh
# report.sh -- regenerate evidence/REPORT.md from the current evidence files.
# Run from kernel/:  ./scripts/report.sh   (after ./scripts/check.sh)
set -eu
cd "$(dirname "$0")/.."
EVID=evidence
BUILD=build
OUT=$EVID/REPORT.md

{
    echo "# Random-source unification -- acceptance evidence"
    echo
    echo "Generated $(date -u +%Y-%m-%dT%H:%M:%SZ) from \`kernel/evidence/\`."
    echo "Regenerate with \`make -C kernel test && kernel/scripts/report.sh\`."
    echo
    echo "## 1. Public API frozen"
    echo
    echo "Header hash (baseline vs now):"
    echo
    echo '```'
    echo "baseline: $(cat baseline/public_api.sha256)"
    echo "now:      $(sha256sum include/rt_kernel.h)"
    echo '```'
    echo
    echo "Compile-time signature assertions:"
    echo
    echo '```'
    cat "$EVID/sigcheck.txt"
    echo '```'
    echo
    echo "Exported symbols (identical between the OFF build and the unified build):"
    echo
    echo '```'
    cat "$EVID/symbols.txt"
    echo '```'
    echo
    echo "## 2. Existing modules untouched"
    echo
    echo "Self-test stdout byte-compared against the frozen baseline:"
    echo
    echo '```'
    for m in bvh sample material; do
        printf '%-9s baseline %s\n' "$m" "$(sha256sum baseline/module_tests/$m.txt | awk '{print $1}')"
        printf '%-9s now      %s\n' "$m" "$(sha256sum $EVID/module_$m.txt | awk '{print $1}')"
    done
    echo '```'
    echo
    echo "Tracked-file diff over bvh/, sample/, material/ (empty means untouched):"
    echo
    echo '```'
    cat "$EVID/module_git_diff.txt"
    echo '```'
    echo
    echo "## 3. Random source property tests"
    echo
    echo '```'
    cat "$EVID/rng_selftest.txt"
    echo '```'
    echo
    echo "## 4. OFF switch: byte-identical to the frozen pre-change output"
    echo
    echo "The OFF build uses the pre-unification per-site generators. For each site,"
    echo "its 256-item stdout was compared with \`cmp\` against the baseline captured"
    echo "before unification:"
    echo
    echo '```'
    for site in 0 1 2; do
        printf 'site%s baseline %s\n' "$site" "$(sha256sum baseline/legacy/site$site.threads1.txt | awk '{print $1}')"
        printf 'site%s current  %s\n' "$site" "$(sha256sum $EVID/legacy_site$site.txt | awk '{print $1}')"
    done
    echo '```'
    echo
    echo "## 5. Thread and schedule independence"
    echo
    echo "Scene fingerprint, all three sites folded into one FNV-1a hash, across"
    echo "thread counts and visit orders:"
    echo
    echo '```'
    cat "$EVID/scene_hash.txt"
    echo '```'
    echo
    echo "Per-site hashes at each thread count are compared in \`check.sh\` section 5."
    echo
    echo "## 6. Error before vs after"
    echo
    echo "Each site estimates an analytic target; 20000 items per run. The verdict"
    echo "requires the mean to be within 5 standard errors of the target."
    echo
    echo '```'
    cat "$EVID/accuracy.txt"
    echo '```'
    echo
    echo "## 7. Mutable generator state"
    echo
    echo '```'
    cat "$EVID/global_state.txt"
    echo '```'
    echo
    echo "## 7b. Negative controls (the checks are not vacuous)"
    echo
    echo "Each control injects one defect into a throwaway copy and requires the"
    echo "corresponding check to fail.  A suite that passes mutants too would"
    echo "measure nothing; these must all report the defect as detected."
    echo
    echo '```'
    cat "$BUILD/neg.out"
    echo '```'
    echo
    echo "## 8. Cross-compiler / cross-optimisation"
    echo
    echo '```'
    cat "$EVID/portability.txt"
    echo '```'
} > "$OUT"

printf 'wrote %s\n' "$OUT"
