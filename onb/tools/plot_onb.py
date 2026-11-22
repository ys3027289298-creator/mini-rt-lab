#!/usr/bin/env python3
"""Plot ONB condition-number curves from sweep.csv / flip.csv.

Usage: plot_onb.py sweep.csv flip.csv out.png
Four panels:
  1. double-precision twist-aligned error vs kappa, with analytic bound;
  2. single-precision same;
  3. raw per-vector error (O(1) geometric twist, kappa-independent);
  4. tangent sign flip across the construction pole.
"""
import glob
import math
import sys
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Validated dataviz categorical palette (adjacent-pair checks all PASS).
BLUE, ORANGE, AQUA, RED = "#2a78d6", "#1baf7a", "#eb6834", "#e34948"
GRID = "#d9d8d3"
INK = "#1d1c1a"
MUTED = "#6b6a66"

# Register a CJK font explicitly: distros ship WenQuanYi as a .ttc, which the
# default font scan may not index.
for pattern in ("/usr/share/fonts/**/wqy-zenhei.ttc",
                "/home/*/.local/share/fonts/wqy-zenhei.ttc",
                "/usr/share/fonts/**/NotoSansCJK*.ttc",
                "/usr/share/fonts/**/SourceHanSans*.otf"):
    hits = glob.glob(pattern, recursive=True)
    if hits:
        font_manager.fontManager.addfont(hits[0])
for name in ("WenQuanYi Zen Hei", "Noto Sans CJK SC", "Noto Sans CJK JP",
             "Source Han Sans SC"):
    if any(name in f.name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.sans-serif"] = [name]
        break
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["text.color"] = INK
plt.rcParams["axes.labelcolor"] = INK
plt.rcParams["xtick.color"] = MUTED
plt.rcParams["ytick.color"] = MUTED


def load(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def fnum(row, key):
    v = float(row[key])
    if v != v or v == 0.0:
        return None          # NaN or exact zero: leave a gap on the log axis
    return v


def ascii_log_label(x, pos=None):
    """Plain '10^-17' labels: the available CJK font has no U+2212 glyph and no
    mathtext fallback is configured, so avoid both unicode minus and mathtext."""
    if x <= 0:
        return ""
    exp = int(round(math.log10(x)))
    return f"10^{exp}" if 10 ** exp == x else ""


def style_axes(ax):
    ax.grid(True, color=GRID, linewidth=0.7, alpha=0.8)
    ax.set_axisbelow(True)


def aligned_panel(ax, rows, prec, title, note_y):
    sub = [r for r in rows if r["precision"] == prec]
    k = [float(r["kappa"]) for r in sub]
    ax.plot(k, [fnum(r, "fixed_aligned") for r in sub], color=ORANGE,
            linestyle="--", linewidth=1.3, marker="o", markersize=3.0,
            label="固定轴构造 · 扭转对齐误差（随 κ 线性放大）")
    ax.plot(k, [fnum(r, "bound") for r in sub], color=BLUE,
            linewidth=1.3, linestyle=":", label="解析上界 C·ε·κ")
    # The stable build has its own kappa, always <= sqrt(2): plot its error
    # against that kappa as a flat band near the left edge, not against the
    # fixed construction's kappa (which would fake a descending curve).
    sk = [float(r["stable_kappa"]) for r in sub]
    se = [fnum(r, "stable_aligned") for r in sub]
    pairs = [(a, b) for a, b in zip(sk, se) if b is not None]
    if pairs:
        ax.plot([p[0] for p in pairs], [p[1] for p in pairs], color=AQUA,
                linewidth=1.3, marker="^", markersize=3.4,
                label="稳定 ONB · 扭转对齐误差（自身 κ ≤ √2，误差恒定有界）")
    ax.axvline(1e4, color=MUTED, linewidth=1.1, linestyle="-.")
    ax.text(1.25e4, note_y, "κ 预算 1e4：超限\n即报告并切换备用构造",
            color=MUTED, fontsize=7.8, va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ascii_log_label))
    ax.set_xlabel("固定构造条件数 κ（绿线使用其自身 κ，落在轴左端）")
    ax.set_ylabel("扭转对齐后的逐向量误差")
    ax.set_title(title, fontsize=10.5)
    style_axes(ax)
    ax.legend(fontsize=7.0, frameon=False, loc="upper left")


def main():
    sweep_path, flip_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    rows = load(sweep_path)

    fig, axes = plt.subplots(1, 4, figsize=(22, 5.3), dpi=150)
    ax1, ax2, ax3, ax4 = axes

    # ---- panels 1-2: twist-aligned numerical error + analytic bound --------
    aligned_panel(ax1, rows, "double",
                  "双精度：误差受 C·ε·κ 上界约束（C=2）", 3e-17)
    aligned_panel(ax2, rows, "float",
                  "单精度：同一条律，ε 放大约 1e9 倍（C=4）", 2e-8)

    # ---- panel 3: raw error = O(1) geometric twist, kappa-independent ------
    for prec, marker in (("double", "o"), ("float", "^")):
        sub = [r for r in rows if r["precision"] == prec]
        k = [float(r["kappa"]) for r in sub]
        ax3.plot(k, [fnum(r, "fixed_raw") for r in sub], color=RED,
                 linewidth=1.4, marker=marker, markersize=3.0,
                 label=f"固定轴 vs 四元数（{prec}）")
        ax3.plot(k, [fnum(r, "stable_raw") for r in sub], color=AQUA,
                 linewidth=1.4, marker=marker, markersize=3.0, linestyle="--",
                 label=f"稳定 ONB vs 四元数（{prec}）")
    ax3.axhline(2.0, color=MUTED, linewidth=1.0, linestyle=":")
    ax3.text(1.4, 2.05, "单位向量逐分量差的平凡上界 2", color=MUTED, fontsize=8.0)
    ax3.set_xscale("log")
    ax3.set_ylim(0, 2.3)
    ax3.set_xlabel("固定构造条件数 κ")
    ax3.set_ylabel("未对齐的原始逐向量误差 ||v_F - v_Q||")
    ax3.set_title("原始误差是 O(1) 几何扭转，不随 κ 发散", fontsize=10.5)
    style_axes(ax3)
    ax3.legend(fontsize=7.4, frameon=False, loc="center right")

    # ---- panel 4: tangent sign flip across the pole ------------------------
    flip = load(flip_path)
    ang_all = [float(r["angle_deg"]) for r in flip]
    fty_all = [fnum(r, "fixed_ty") for r in flip]
    sty_all = [float(r["stable_ty"]) for r in flip]
    fk_all = [float(r["fixed_kappa"]) for r in flip]
    sk_all = [float(r["stable_kappa"]) for r in flip]
    mask = [abs(a) <= 1.0 for a in ang_all]
    ang = [a for a, ok in zip(ang_all, mask) if ok]
    fty = [v for v, ok in zip(fty_all, mask) if ok]
    sty = [v for v, ok in zip(sty_all, mask) if ok]
    fk = [v for v, ok in zip(fk_all, mask) if ok]
    sk = [v for v, ok in zip(sk_all, mask) if ok]

    ax4.plot(ang, fty, color=RED, linewidth=1.8,
             label="固定轴构造：切向量 t_y（跨极轴变号）")
    ax4.plot(ang, sty, color=AQUA, linewidth=1.8,
             label="稳定 ONB：切向量 t_y（保持连续）")
    ax4.axvline(0.0, color=MUTED, linewidth=1.0, linestyle="--")
    ax4.scatter([0], [1], color=AQUA, s=28, zorder=5)
    ax4.set_xlabel("法线在 xz 平面内的倾角（度），扫过 +z 极轴")
    ax4.set_ylabel("切向量 t 的 y 分量")
    ax4.set_title("极轴附近切向翻转——余弦采样方向被镜像的根因", fontsize=10.5)
    style_axes(ax4)
    fk_finite = max(x for x in fk if x == x)
    ax4.annotate(f"固定构造在 0° 未定义（邻点 κ={fk_finite:.1e}）\n"
                 f"稳定构造 κ 恒为 {max(sk):.4f} ≤ √2，极轴处仍给出框架",
                 xy=(0, -1), xytext=(0.12, -0.68),
                 fontsize=7.8, color=MUTED,
                 arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0))
    ax4.legend(fontsize=7.8, frameon=False, loc="center right")
    ax4.set_ylim(-1.28, 1.28)

    fig.tight_layout(w_pad=2.4)
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    print(f"plot_written={out_path}")


if __name__ == "__main__":
    main()
