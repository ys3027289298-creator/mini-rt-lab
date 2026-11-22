"""Render the volume and residual curves from records.json.

Two single-axis panels (never a dual axis):
  A. residual vs rank for the block codec, with its Eckart-Young lower bound
     drawn underneath (they coincide: a rank-r SVD leaves exactly the discarded
     energy), plus the unblocked global SVD as a comparison series;
  B. storage volume vs rank, with the lossless reference as a horizontal rule.

Palette: validated categorical slots 1-2 (blue, orange); theory/reference ink is
muted chrome, not a series.
"""
import json
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FuncFormatter

# --- palette (light surface), validated with dataviz/scripts/validate_palette.js
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"     # slot 1: block codec
ORANGE = "#eb6834"   # slot 2: global codec

_FONT_CANDIDATES = [
    "/home/ys/.local/share/fonts/wqy-zenhei.ttc",
    "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
]


def _setup_font():
    for path in _FONT_CANDIDATES:
        if os.path.exists(path):
            font_manager.fontManager.addfont(path)
            name = font_manager.FontProperties(fname=path).get_name()
            plt.rcParams["font.family"] = [name, "DejaVu Sans"]
            break
    else:
        plt.rcParams["font.family"] = ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False


def _log_tick(y, _pos):
    # plain-text scientific label; avoids mathtext so the CJK font never needs
    # to supply U+2212 (the minus mathtext would emit)
    if y <= 0:
        return ""
    exp = int(np.floor(np.log10(y)))
    mant = y / 10.0 ** exp
    if abs(mant - 1.0) < 1e-9:
        return f"1e{exp}"
    return f"{mant:g}e{exp}" if abs(exp) < 3 else f"{mant:g}e{exp}"


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=9, length=3, width=0.8)
    ax.grid(True, color=GRID, linewidth=0.8, linestyle="-", alpha=0.9)
    ax.set_axisbelow(True)


def main(outdir):
    _setup_font()
    with open(os.path.join(outdir, "records.json")) as f:
        rec = json.load(f)

    block = [r for r in rec["block"]["records"] if r["rank"] > 0]
    glob = rec["global"]["records"]
    lossless = rec["lossless"]
    H, W = rec["frame"]["height"], rec["frame"]["width"]
    N = H * W * rec["frame"]["channels"]

    rb = np.array([r["rank"] for r in block])
    rmse_b = np.array([r["rmse"] for r in block])
    bound_b = np.sqrt(np.array([r["bound_sq"] for r in block]) / N)
    vol_b = np.array([r["bytes"] for r in block]) / 1024.0

    rg = np.array([r["rank"] for r in glob])
    rmse_g = np.array([r["rmse"] for r in glob])
    vol_g = np.array([r["bytes"] for r in glob]) / 1024.0

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11.6, 4.5), facecolor=SURFACE, constrained_layout=True
    )
    _style(ax1)
    _style(ax2)

    # ---- panel A: residual vs rank, bound drawn underneath as a translucent floor
    pos = bound_b > 0
    ax1.plot(rb[pos], bound_b[pos], color=MUTED, linewidth=4.5, alpha=0.55,
             solid_capstyle="round", label="Eckart–Young 下界（块 SVD）", zorder=1)
    ax1.plot(rb, rmse_b, color=BLUE, linewidth=1.8, marker="o", markersize=5,
             label="块 SVD 实测残差", zorder=3)
    ax1.plot(rg, rmse_g, color=ORANGE, linewidth=1.8, marker="s", markersize=4.6,
             linestyle=(0, (5, 2)), label="全局 SVD 实测残差", zorder=2)
    ax1.set_yscale("log")
    ax1.yaxis.set_major_formatter(FuncFormatter(_log_tick))
    ax1.yaxis.set_minor_formatter(FuncFormatter(lambda *_: ""))
    ax1.set_xlabel("秩 r", color=INK2, fontsize=10)
    ax1.set_ylabel("RMSE（线性 HDR，按像素归一）", color=INK2, fontsize=10)
    ax1.set_title("残差曲线：实测残差贴着理论下界，从不低于它", color=INK,
                  fontsize=11.5, loc="left", pad=8)
    ratio = rmse_b[pos] / bound_b[pos]
    ax1.annotate(
        f"实测/下界 ∈ [{ratio.min():.3f}, {ratio.max():.3f}]\n"
        "灰色为下界，蓝色落在其上 → 秩 r 截断即最优",
        (rb[len(rb) // 3], rmse_b[len(rb) // 3]),
        textcoords="offset points", xytext=(26, -46), color=INK2,
        fontsize=8.5, ha="left",
        arrowprops=dict(arrowstyle="-", color=AXIS, linewidth=0.9))
    ax1.annotate(f"r={int(rb[0])}: {rmse_b[0]:.3g}", (rb[0], rmse_b[0]),
                 textcoords="offset points", xytext=(6, -12), color=BLUE, fontsize=8.5)
    ax1.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower left")

    # ---- panel B: volume vs rank
    ax2.plot(rb, vol_b, color=BLUE, linewidth=1.8, marker="o", markersize=5,
             label="块 SVD", zorder=3)
    ax2.plot(rg, vol_g, color=ORANGE, linewidth=1.8, marker="s", markersize=4.6,
             linestyle=(0, (5, 2)), label="全局 SVD", zorder=3)
    lossless_kb = lossless["bytes"] / 1024.0
    ax2.axhline(lossless_kb, color=MUTED, linewidth=1.4, zorder=2)
    ax2.annotate(f"无损参照 {lossless_kb:.1f} KiB", (rb[0], lossless_kb),
                 textcoords="offset points", xytext=(2, 5), color=INK2, fontsize=8.5)
    ax2.set_xlabel("秩 r", color=INK2, fontsize=10)
    ax2.set_ylabel("体积 KiB", color=INK2, fontsize=10)
    ax2.set_title("体积曲线：低秩显著小于无损，高秩反超", color=INK,
                  fontsize=11.5, loc="left", pad=8)
    ax2.annotate(f"{vol_b[0]:.1f} KiB", (rb[0], vol_b[0]),
                 textcoords="offset points", xytext=(10, 8), color=BLUE, fontsize=8.5)
    ax2.annotate(f"{vol_b[-1]:.1f} KiB", (rb[-1], vol_b[-1]),
                 textcoords="offset points", xytext=(-6, 8), color=BLUE, fontsize=8.5,
                 ha="right")
    ax2.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper left")

    fig.suptitle(
        f"块大小 {rec['block']['block_size']} · 因子 {rec['block']['factor_bits']} bit · "
        f"{rec['frame']['width']}×{rec['frame']['height']} · "
        f"下界校验：{'通过' if rec['bounds_ok'] else '失败'}",
        color=INK2, fontsize=9.5, x=0.012, ha="left",
    )
    path = os.path.join(outdir, "curves.png")
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "out")
