"""Compression experiment driver.

Renders the target frame, sweeps ranks for the block codec, and writes:
  out/frame.npz        - the source frame (float32, HDR linear)
  out/records.json     - every measurement (residual, bound, volume, flags)
  out/report.md        - tables for the volume and residual/bound curves
  out/curves.png       - the two curves (volume, residual vs. Eckart-Young bound)

Exit status is non-zero if the measured residual falls below the theoretical
lower bound, i.e. if the implementation is biased.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "src"))

import rt
import lowrank


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=240)
    ap.add_argument("--spp", type=int, default=64)
    ap.add_argument("--block", type=int, default=16)
    ap.add_argument("--bits", type=int, default=16)
    ap.add_argument("--rank-max", type=int, default=16)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--outdir", default=os.path.join(_HERE, "..", "out"))
    a = ap.parse_args(argv)

    if not (4 <= a.block <= 64):
        ap.error("--block must be in [4, 64]")
    if a.width % a.block or a.height % a.block:
        ap.error(f"{a.width}x{a.height} not divisible by block {a.block}")
    if a.bits not in (16, 32):
        ap.error("--bits must be 16 or 32")

    outdir = os.path.abspath(a.outdir)
    os.makedirs(outdir, exist_ok=True)
    energy = None

    print(f"render {a.width}x{a.height} spp={a.spp} ...", flush=True)
    img = rt.render(a.width, a.height, spp=a.spp, seed=a.seed).astype(np.float32)
    vol = img.astype(np.float64)
    energy = lowrank.volume_energy(vol)
    print(f"  mean={vol.mean():.4f} max={vol.max():.4f} "
          f"min={vol.min():.4f}", flush=True)

    # lossless reference: also the destination of the rank-0 fallback
    _, lossless = lowrank.compress_lossless(vol)
    lossless = lowrank.finish_stats(vol, lossless, vol)
    lossless["energy_sq"] = energy
    lossless_bits = lossless["bytes"]

    chans = lowrank.prepare_channels(vol, a.block)
    rmax = min(a.rank_max, a.block)
    block_recs = []
    for r in range(0, rmax + 1):
        rec, st = lowrank.compress_blocks(vol, chans, r, a.bits)
        st["lossless_bytes"] = lossless_bits
        st["energy_sq"] = energy
        lowrank.finish_stats(rec, st, vol)
        block_recs.append(st)

    # unblocked SVD at a few ranks, for comparison only
    global_recs = []
    for r in sorted({1, min(4, min(a.height, a.width)),
                     min(16, min(a.height, a.width))}):
        rec, st = lowrank.compress_global(vol, r)
        st["lossless_bytes"] = lossless_bits
        st["energy_sq"] = energy
        lowrank.finish_stats(rec, st, vol)
        global_recs.append(st)

    checked = [r for r in block_recs if r["rank"] > 0] + global_recs
    ok, worst = lowrank.check_bounds(checked)
    ratios = [r["squared_error"] / r["bound_sq"]
              for r in checked if r["bound_sq"] > 0]
    ratio_lo = min(ratios) if ratios else 1.0
    ratio_hi = max(ratios) if ratios else 1.0
    # worst margin as a relative quantity is the legible form of "how far above
    # the bound" -- an absolute squared-error residual is meaningless on its own
    rel_margin = worst / max(energy, 1e-30)

    records = {
        "frame": {
            "width": a.width, "height": a.height, "spp": a.spp,
            "seed": a.seed, "channels": 3, "dtype": "float32",
            "source_format": "HDR linear RGB (no tone map)",
        },
        "lossless": lossless,
        "block": {"block_size": a.block, "factor_bits": a.bits,
                  "records": block_recs},
        "global": {"records": global_recs},
        "bounds_ok": ok,
        "worst_margin": worst,
        "bound_ratio_min": ratio_lo,
        "bound_ratio_max": ratio_hi,
        "worst_margin_rel": rel_margin,
    }
    with open(os.path.join(outdir, "records.json"), "w") as f:
        json.dump(records, f, indent=2)
    np.savez_compressed(os.path.join(outdir, "frame.npz"), image=img)

    _write_report(outdir, records)
    print(f"bounds_ok={int(ok)} worst_margin_rel={rel_margin:.3g} "
          f"measured/bound=[{ratio_lo:.6f},{ratio_hi:.6f}]")
    print(f"lossless {lossless['bytes']:,} B; "
          f"block r={rmax} {block_recs[-1]['bytes']:,} B")

    try:
        sys.path.insert(0, os.path.join(_HERE, "..", "tools"))
        import make_chart
        make_chart.main(outdir)
    except Exception as e:  # chart is reporting only, never fatal
        print(f"[chart skipped: {e}]", file=sys.stderr)

    if not ok:
        print("FAIL: measured residual below the Eckart-Young lower bound; "
              "the implementation is biased", file=sys.stderr)
        return 2
    return 0


def _fmt(v, nd=3):
    if isinstance(v, float) and (v != v or v in (float("inf"), float("-inf"))):
        return "∞" if v > 0 else str(v)
    return f"{v:.{nd}f}" if isinstance(v, float) else str(v)


def _write_report(outdir, rec):
    L = []
    L.append("# 低秩压缩实验报告\n")
    f = rec["frame"]
    L.append(f"- 帧: {f['width']}x{f['height']} RGB，spp={f['spp']}，seed={f['seed']}，"
             f"线性 HDR（{f['source_format']}）")
    L.append(f"- 无损参照体积: **{rec['lossless']['bytes']:,} B**"
             f"（未压缩 {rec['lossless']['raw_bytes']:,} B）")
    L.append(f"- 块大小 {rec['block']['block_size']}，因子位深 "
             f"{rec['block']['factor_bits']} bit/分量")
    L.append(f"- 理论下界校验: **{'通过' if rec['bounds_ok'] else '失败'}**；"
             f"各秩实测/下界 ∈ [{rec['bound_ratio_min']:.6f}, "
             f"{rec['bound_ratio_max']:.6f}]，最小相对余量 "
             f"{rec['worst_margin_rel']:.3g}（相对总能量）\n")

    L.append("## 残差与理论下界\n")
    L.append("Eckart-Young：任意秩 r 近似的平方误差 ≥ 该块被丢弃奇异值平方和。"
             "`实测/下界` ≥ 1 表示实测残差高于理论下界（因含量化误差），自洽。\n")
    L.append("| 模式 | 秩 r | RMSE | PSNR dB | 实测平方误差 | 理论下界 | 实测/下界 |")
    L.append("|---|---|---|---|---|---|---|")
    for r in rec["block"]["records"] + rec["global"]["records"]:
        if r["mode"] == "block" and r["rank"] == 0:
            continue
        ratio = (r["squared_error"] / r["bound_sq"]) if r["bound_sq"] > 0 else float("inf")
        L.append(f"| {r['mode']} | {r['rank']} | {_fmt(r['rmse'])} | "
                 f"{_fmt(r['psnr_db'], 1)} | {r['squared_error']:.6g} | "
                 f"{r['bound_sq']:.6g} | {_fmt(ratio)} |")

    L.append("\n## 体积\n")
    L.append("体积为真实无损编码后的字节数（比特平面 + deflate），逐块校验往返一致。\n")
    L.append("| 模式 | 秩 r | 体积 B | 相对无损 | 相对原始 | 每像素 B |")
    L.append("|---|---|---|---|---|---|")
    for r in rec["block"]["records"] + rec["global"]["records"]:
        L.append(f"| {r['mode']} | {r['rank']} | {r['bytes']:,} | "
                 f"{_fmt(r['ratio_vs_lossless'], 3)}x | {_fmt(r['ratio_vs_raw'], 3)}x | "
                 f"{_fmt(r['bytes_per_pixel'], 4)} |")

    L.append("\n## 退化：秩 0\n")
    r0 = [x for x in rec["block"]["records"] if x["rank"] == 0][0]
    L.append(f"秩 0 不携带任何图像信息，自动退回无损路径：体积 {r0['bytes']:,} B，"
             f"平方误差 {r0['squared_error']:.6g}，fallback={r0['fallback']}，"
             f"原因 `{r0.get('fallback_reason', '')}`。")

    with open(os.path.join(outdir, "report.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
