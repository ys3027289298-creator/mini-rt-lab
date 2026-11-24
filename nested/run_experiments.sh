#!/usr/bin/env bash
# Full nested-sphere experiment.
#
# For each inner-sphere albedo rho:
#   1. derive penetration cap K and energy threshold tau from rho
#      (K = smallest integer with rho^K <= 1/255, tau = rho^K);
#   2. render the legacy one-query tracer (buggy) and the continuous walker
#      with the derived limits (walk, the reference answer);
#   3. SPP=1 pair so hard chunks show without MC averaging;
#   4. cap sweeps under three threshold regimes:
#        derived tau (correct calibration),
#        tau too HIGH  = rho^2 (cuts paths while energy is still visible;
#                                raising the cap cannot help),
#        tau too LOW   = LSB/100 (never fires at visible energy; the cap
#                                  alone binds, at the cost of more work).
# Every rho group is analyzed against its walk reference.
set -euo pipefail

cd "$(dirname "$0")"
W=320; H=320; SPP=64; SEED=20261124
OUT=renders/exp
mkdir -p "$OUT"
LOG="$OUT/run.log"
: > "$LOG"

run() {
  echo "+ nested $*" | tee -a "$LOG"
  ./nested "$@" 2>&1 | tee -a "$LOG"
}

for RHO in 0.25 0.75; do
  R="$OUT/rho$RHO"
  mkdir -p "$R"

  # Derive K/tau from rho via the binary's own derivation (-1 -1).
  LINE=$(./nested "$W" "$H" "$SPP" "$RHO" walk -1 -1 "$SEED" "$R/probe.ppm" | head -1)
  rm -f "$R/probe.ppm" "$R/probe.mask.ppm"
  K=$(echo "$LINE" | sed -n 's/.*K=\([0-9]*\).*/\1/p')
  TAU=$(echo "$LINE" | sed -n 's/.*tau=\([0-9.]*\).*/\1/p')
  # Fires immediately after the first penetration, while surviving energy
  # rho (0.25 or 0.75) is far above the LSB: a visibly wrong calibration.
  TAU_HIGH=$(python3 -I -c "print('%.8f' % ($RHO * 0.5))")
  TAU_LOW=$(python3 -I -c "print('%.8f' % (1.0/255.0/100.0))")
  echo "rho=$RHO derived K=$K tau=$TAU | high=$TAU_HIGH low=$TAU_LOW" | tee -a "$LOG"

  # Headline renders.
  run "$W" "$H" "$SPP" "$RHO" buggy "$K" "$TAU" "$SEED" "$R/buggy.ppm"
  run "$W" "$H" "$SPP" "$RHO" walk  "$K" "$TAU" "$SEED" "$R/walk.ppm"

  # SPP=1 pair: hard chunks pixel for pixel, plus SPP=1 threshold sweeps so
  # threshold-mis-set effects show as actual black pixels instead of a
  # sub-code-value average bias.
  run "$W" "$H" 1 "$RHO" buggy "$K" "$TAU" "$SEED" "$R/buggy_spp1.ppm"
  run "$W" "$H" 1 "$RHO" walk  "$K" "$TAU" "$SEED" "$R/walk_spp1.ppm"
  for CK in 1 2 4 8 16; do
    run "$W" "$H" 1 "$RHO" walk "$CK" "$TAU_HIGH" "$SEED" \
        "$R/walk_spp1_tauHigh_k$CK.ppm"
    run "$W" "$H" 1 "$RHO" walk "$CK" "$TAU_LOW" "$SEED" \
        "$R/walk_spp1_tauLow_k$CK.ppm"
  done

  for CK in 1 2 4 8 16; do
    run "$W" "$H" "$SPP" "$RHO" walk "$CK" "$TAU"      "$SEED" \
        "$R/walk_tauDerived_k$CK.ppm"
    run "$W" "$H" "$SPP" "$RHO" walk "$CK" "$TAU_HIGH" "$SEED" \
        "$R/walk_tauHigh_k$CK.ppm"
    run "$W" "$H" "$SPP" "$RHO" walk "$CK" "$TAU_LOW"  "$SEED" \
        "$R/walk_tauLow_k$CK.ppm"
  done

  # SPP=64 group (reference: the derived-limit walk render).
  TESTS_64=$(ls "$R" | grep '\.ppm$' | grep -v '\.mask\.ppm' \
             | grep -v '^walk.ppm' | grep -v '_spp1' | sed "s|^|$R/|")
  python3 -I tools/analyze.py "$R/analysis64" "$R/walk.ppm" $TESTS_64

  # SPP=1 group with its own SPP=1 reference.
  TESTS_1=$(ls "$R" | grep '_spp1.*\.ppm$' | grep -v '\.mask\.ppm' \
            | grep -v '^walk_spp1.ppm' | sed "s|^|$R/|")
  python3 -I tools/analyze.py "$R/analysis1" "$R/walk_spp1.ppm" $TESTS_1
done
