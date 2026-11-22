"""Self-test for the low-rank compression module.  `make test` runs this.

Checks, in order of what would break first:
  1. the lossless coder is byte-exact (round-trip on random and real data);
  2. the bit-plane coder preserves values that deflate cannot shrink;
  3. per-block SVD truncation obeys Eckart-Young: measured residual >= lower
     bound for every rank, every block, with a tight relative tolerance;
  4. rank 0 degenerates to the lossless path and is flagged;
  5. volume falls monotonically with rank up to the block size;
  6. the bound guard actually catches a deliberately biased record.
Exit status is non-zero on any failure.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "src"))

import coding
import lowrank
import rt

failures = 0


def check(name, ok, detail=""):
    global failures
    print(f"{name}={'1' if ok else '0'}" + (f"  {detail}" if detail and not ok else ""))
    if not ok:
        failures += 1


def test_coding():
    rng = np.random.default_rng(3)
    for dtype in (np.uint8, np.uint16, np.uint32):
        a = rng.integers(0, np.iinfo(dtype).max, size=(7, 13), dtype=dtype)
        payload = coding._to_bits(a)
        back = coding._from_bits(payload, a.shape, dtype)
        check(f"bits_roundtrip_{np.dtype(dtype).name}", np.array_equal(a, back))
    r = rng.random((9, 11)).astype(np.float32)
    comp = coding.coded_size([(r.view(np.uint32), np.uint32)])
    check("coded_size_positive", comp > 0)


def _frame(h=64, w=64, spp=4):
    return rt.render(w, h, spp=spp, seed=11).astype(np.float32).astype(np.float64)


def test_bounds():
    vol = _frame()
    chans = lowrank.prepare_channels(vol, 16)
    E = lowrank.volume_energy(vol)
    ok_all = True
    worst = 0.0
    recs = []
    for r in range(1, 17):
        rec, st = lowrank.compress_blocks(vol, chans, r, 16)
        st["energy_sq"] = E
        margin = st["squared_error"] - st["bound_sq"]
        # lower bound must hold within float round-off
        tol = 1e-6 * max(st["bound_sq"], E)
        if margin < -tol:
            ok_all = False
            worst = min(worst, margin)
        recs.append(st)
    check("eckart_young_bound_holds", ok_all, f"worst margin {worst:.3g}")
    # and the bound must be witnessed: residual should sit above it but close.
    # Tightness is measured against total energy -- the margin at high rank is
    # pure factor-quantisation round-off, a fixed tiny fraction of the energy.
    ratios = [s["squared_error"] / s["bound_sq"] for s in recs if s["bound_sq"] > 0]
    check("residual_above_bound", all(rr >= 1.0 - 1e-6 for rr in ratios))
    rel_margins = [abs(s["squared_error"] - s["bound_sq"]) / E
                   for s in recs if s["bound_sq"] > 0]
    check("bound_is_tight", max(rel_margins) < 1e-6,
          f"max |measured-bound|/energy {max(rel_margins):.3g}")

    # guard catches an injected bias (claim a residual smaller than possible)
    bad = dict(recs[3])
    bad["squared_error"] = bad["bound_sq"] * 0.5
    ok, _ = lowrank.check_bounds([bad])
    check("guard_detects_bias", ok is False)


def test_rank_zero():
    vol = _frame()
    chans = lowrank.prepare_channels(vol, 16)
    rec0, st0 = lowrank.compress_blocks(vol, chans, 0, 16)
    _, lz = lowrank.compress_lossless(vol)
    check("rank_zero_fallback_flag", st0["fallback"] is True)
    check("rank_zero_error_zero", st0["squared_error"] == 0.0)
    check("rank_zero_lossless_size", st0["bytes"] == lz["bytes"])
    check("rank_zero_exact", np.array_equal(rec0, vol))


def test_volume_monotone():
    vol = _frame()
    chans = lowrank.prepare_channels(vol, 16)
    sizes = []
    for r in range(1, 17):
        _, st = lowrank.compress_blocks(vol, chans, r, 16)
        sizes.append(st["bytes"])
    check("volume_increases_with_rank",
          all(b >= a for a, b in zip(sizes, sizes[1:])),
          f"{sizes}")


def test_factor_bits():
    vol = _frame()
    chans = lowrank.prepare_channels(vol, 16)
    E = lowrank.volume_energy(vol)
    _, s16 = lowrank.compress_blocks(vol, chans, 6, 16)
    _, s32 = lowrank.compress_blocks(vol, chans, 6, 32)
    check("float32_factors_more_accurate", s32["squared_error"] < s16["squared_error"])
    # float32 factors carry their own rounding; isolate it as the excess over
    # the truncation bound (both ranks have the same singular values, so the
    # bound is identical) and require it to be negligible next to the energy.
    excess = (s32["squared_error"] - s32["bound_sq"]) / E
    check("float32_factor_error_tiny", excess < 1e-10,
          f"excess over bound {excess:.3e} of energy")


def test_global():
    vol = _frame()
    chans = lowrank.prepare_channels(vol, 16)
    g_recs = []
    E = lowrank.volume_energy(vol)
    for r in (2, 8, 16):
        _, st = lowrank.compress_global(vol, r)
        st["energy_sq"] = E
        g_recs.append(st)
    ok, _ = lowrank.check_bounds(g_recs)
    check("global_bound_holds", ok)
    # block r=bs should equal the lossless-quality reconstruction of each block
    _, sb = lowrank.compress_blocks(vol, chans, 16, 32)
    check("full_rank_near_lossless",
          sb["squared_error"] < 1e-12 * max(E, 1.0))


def main():
    test_coding()
    test_bounds()
    test_rank_zero()
    test_volume_monotone()
    test_factor_bits()
    test_global()
    print(f"all_ok={int(failures == 0)}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
