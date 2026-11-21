#!/usr/bin/env python3
"""Plot the fixed-budget resolution/spp trade-off from the C module's CSV.

The C binary is the single source of truth: it enumerates feasible combos,
the exhaustive optimum, the fast pick, and the continuous relaxation. This
script only parses its `csv` output and draws:

  budget_curves.png   (a) error vs spp, one line per resolution, plus the
                          continuous lower bound E*_c and opt/pick markers
                      (b) per-resolution best error vs pixel count
  budget_bounds.png   across a budget sweep, two panels:
                      top    absolute pick-vs-exhaustive gap vs its certified
                             instance upper bound (log scale)
                      bottom relative gap / bound in % of E*_c (linear, so the
                             zero gaps where the pick hits the optimum show)

Usage: python3 plot_budget.py [B] [a] [v] [outdir]
"""
import csv
import math
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- validated reference palette (dataviz skill, light mode) ---
SURFACE = "#fcfcfb"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MARKERS = ["o", "s", "^", "D", "v", "P"]
BOUND_HUE = "#e34948"          # slot 8 red, certified bound
LB_HUE = "#6250d6"             # slot 7 violet, continuous lower bound
SHADE = "#e1e0d9"              # infeasible region shade

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 10,
    "text.color": INK,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


def run_csv(binary: Path, B: int, a: float, v: float):
    out = subprocess.run(
        [str(binary), "csv", str(B), str(a), str(v)],
        check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    meta = {}
    for token in out[0].lstrip("# ").split():
        key, val = token.split("=", 1)
        meta[key] = val
    rows = list(csv.DictReader(line for line in out[1:] if line))
    for r in rows:
        for k in ("w", "h", "spp"):
            r[k] = int(r[k])
        for k in ("pixels", "used"):
            r[k] = int(r[k])
        r["err"] = float(r["err"])
    return meta, rows


def plot_curves(meta, rows, B, a, v, path: Path):
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.5, 5.4), gridspec_kw={"width_ratios": [1.55, 1.0]}
    )
    fig.suptitle(
        f"Fixed total budget B = {B:,} samples   (a = {a:g}, v = {v:g}, "
        f"E = a/pixels + v/spp)",
        fontsize=12, color=INK,
    )

    pixels_levels = sorted({r["pixels"] for r in rows})
    color_of = {p: SERIES[i % len(SERIES)] for i, p in enumerate(pixels_levels)}
    mark_of = {p: MARKERS[i % len(MARKERS)] for i, p in enumerate(pixels_levels)}
    exh = next(r for r in rows if r["role"] in ("exhaustive", "both"))
    pick = next(r for r in rows if r["role"] in ("pick", "both"))
    ec = float(meta["ec"])

    # (a) error curves, one per resolution --------------------------------
    for px in pixels_levels:
        pts = sorted((r for r in rows if r["pixels"] == px),
                     key=lambda r: r["spp"])
        c = color_of[px]
        side = int(math.sqrt(px))
        ax1.plot([r["spp"] for r in pts], [r["err"] for r in pts],
                 color=c, linewidth=2, marker=mark_of[px], markersize=6,
                 markerfacecolor=c, markeredgecolor=SURFACE,
                 markeredgewidth=1.5,
                 label=f"{side}x{side} ({px:,} px)")

    ax1.axhline(ec, color=LB_HUE, linewidth=1.6, linestyle=(0, (6, 4)))
    ax1.annotate(f"continuous lower bound E*_c = {ec:.4g}",
                 xy=(0.02, ec), xycoords=("axes fraction", "data"),
                 xytext=(0, 8), textcoords="offset points",
                 ha="left", va="bottom", fontsize=9, color=INK2)

    ax1.scatter([exh["spp"]], [exh["err"]], s=200, facecolors="none",
                edgecolors=INK, linewidths=1.8, zorder=6)
    ax1.annotate("exhaustive optimum\n= fast pick" if exh is pick
                 else "exhaustive optimum",
                 xy=(exh["spp"], exh["err"]), xytext=(-18, -18),
                 textcoords="offset points", ha="right", va="top",
                 fontsize=9, color=INK,
                 arrowprops=dict(arrowstyle="-", color=INK2, lw=0.9))
    if pick is not exh:
        ax1.scatter([pick["spp"]], [pick["err"]], s=160, marker="*",
                    facecolors=BOUND_HUE, edgecolors=SURFACE,
                    linewidths=1.2, zorder=6)
        ax1.annotate("fast pick", xy=(pick["spp"], pick["err"]),
                     xytext=(8, 12), textcoords="offset points",
                     fontsize=9, color=INK,
                     arrowprops=dict(arrowstyle="-", color=INK2, lw=0.9))

    ax1.set_xscale("log", base=2)
    ax1.set_yscale("log")
    ax1.set_xmargin(0.12)
    ax1.set_ymargin(0.22)
    ax1.set_xlabel("samples per pixel (spp)")
    ax1.set_ylabel("error E")
    ax1.set_title("(a) error curves per resolution (feasible combos only)",
                  fontsize=10.5, color=INK)
    ax1.legend(title="resolution", fontsize=8.5, title_fontsize=9,
               ncol=2, loc="upper right")

    # (b) best error attainable at each resolution ------------------------
    best = []
    for px in pixels_levels:
        r = min((r for r in rows if r["pixels"] == px),
                key=lambda r: r["err"])
        best.append(r)
    ax2.plot([r["pixels"] for r in best], [r["err"] for r in best],
             color=SERIES[0], linewidth=2, marker="o", markersize=7,
             markerfacecolor=SERIES[0], markeredgecolor=SURFACE,
             markeredgewidth=1.5)
    ax2.axhline(ec, color=LB_HUE, linewidth=1.6, linestyle=(0, (6, 4)))
    ax2.scatter([exh["pixels"]], [exh["err"]], s=200, facecolors="none",
                edgecolors=INK, linewidths=1.8, zorder=6)
    for i, r in enumerate(best):
        # alternate label side so adjacent annotations never collide
        dy = 9 if i % 2 == 0 else -13
        ax2.annotate(f"spp {r['spp']}", xy=(r["pixels"], r["err"]),
                     xytext=(0, dy), textcoords="offset points",
                     ha="center", fontsize=8, color=MUTED)
    ax2.set_xscale("log", base=2)
    ax2.set_yscale("log")
    ax2.set_xmargin(0.10)
    ax2.set_xlabel("pixels P")
    ax2.set_ylabel("minimum error at this resolution")
    ax2.set_title("(b) per-resolution best vs continuous bound",
                  fontsize=10.5, color=INK)

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_bounds(binary: Path, a: float, v: float, path: Path):
    budgets = [2 ** k for k in range(5, 25)]          # 32 .. 16M; 32 < min cost 64
    feasible_b, gaps, bounds, rel_gaps, rel_bounds, infeasible = [], [], [], [], [], []
    for B in budgets:
        meta, rows = run_csv(binary, B, a, v)
        if meta["feasible"] != "1":
            infeasible.append(B)
            continue
        exh = next(r for r in rows if r["role"] in ("exhaustive", "both"))
        pick = next(r for r in rows if r["role"] in ("pick", "both"))
        ec = float(meta["ec"])
        gap = pick["err"] - exh["err"]
        bound = pick["err"] - ec
        feasible_b.append(B)
        gaps.append(gap)
        bounds.append(bound)
        rel_gaps.append(100.0 * gap / ec)
        rel_bounds.append(100.0 * bound / ec)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(10.5, 7.6), sharex=True,
        gridspec_kw={"height_ratios": [1.25, 1.0]},
    )
    fig.suptitle(
        f"Fast pick vs exhaustive optimum over a budget sweep "
        f"(a = {a:g}, v = {v:g})",
        fontsize=12, color=INK,
    )

    shade_hi = min(feasible_b)
    for ax in (ax_top, ax_bot):
        if infeasible:
            ax.axvspan(min(budgets) / 1.25, shade_hi, color=SHADE,
                       alpha=0.75, zorder=0)
    ax_top.annotate("no feasible solution\n(empty combo set)",
                    xy=(math.sqrt((min(budgets) / 1.25) * shade_hi), 0.97),
                    xycoords=("data", "axes fraction"),
                    ha="center", va="top", fontsize=9, color=INK2)

    # top: absolute deviation, log scale; zero gaps use the axis floor ----
    floor = min(b for b in bounds if b > 0) * 0.35
    gap_plot = [g if g > 0 else floor for g in gaps]
    ax_top.plot(feasible_b, bounds, color=BOUND_HUE, linewidth=2,
                marker="D", markersize=6, markerfacecolor=BOUND_HUE,
                markeredgecolor=SURFACE, markeredgewidth=1.5,
                label="certified upper bound  E_pick − E*_c")
    ax_top.plot(feasible_b, gap_plot, color=SERIES[0], linewidth=2,
                marker="o", markersize=7, markerfacecolor=SERIES[0],
                markeredgecolor=SURFACE, markeredgewidth=1.5,
                label="observed deviation  E_pick − E_exhaustive "
                      "(zeros plotted on floor)")
    ax_top.set_yscale("log")
    ax_top.set_ylabel("absolute error deviation")
    ax_top.set_title("absolute deviation and certified bound",
                     fontsize=10.5, color=INK, loc="left")
    ax_top.legend(loc="upper right", fontsize=9)

    # bottom: relative deviation, linear scale, zeros are real -------------
    ax_bot.plot(feasible_b, rel_bounds, color=BOUND_HUE, linewidth=2,
                marker="D", markersize=6, markerfacecolor=BOUND_HUE,
                markeredgecolor=SURFACE, markeredgewidth=1.5,
                label="bound as % of E*_c")
    ax_bot.plot(feasible_b, rel_gaps, color=SERIES[0], linewidth=2,
                marker="o", markersize=7, markerfacecolor=SERIES[0],
                markeredgecolor=SURFACE, markeredgewidth=1.5,
                label="observed gap as % of E*_c (0 = pick hits optimum)")
    ax_bot.axhline(125.0, color=MUTED, linewidth=1.3,
                   linestyle=(0, (6, 4)))
    ax_bot.annotate("structural guarantee 125%",
                    xy=(feasible_b[-1], 125.0), xytext=(0, 5),
                    textcoords="offset points", ha="right", va="bottom",
                    fontsize=8.5, color=INK2)
    ax_bot.set_ylim(bottom=0)
    ax_bot.set_ylabel("relative deviation (%)")
    ax_bot.set_title("relative deviation (linear scale: zero gaps visible)",
                     fontsize=10.5, color=INK, loc="left")
    ax_bot.legend(loc="upper right", fontsize=9)

    ax_bot.set_xscale("log", base=2)
    ax_bot.set_xmargin(0.06)
    ax_bot.set_xlabel("total budget B (samples)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 1_048_576
    a = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
    v = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    outdir = Path(sys.argv[4]) if len(sys.argv) > 4 else Path("../artifacts")
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    binary = (Path(__file__).resolve().parent.parent / "budget")
    meta, rows = run_csv(binary, B, a, v)
    if meta["feasible"] != "1":
        print(f"no feasible solution for B={B}; skipping curve plot "
              f"(bound sweep still drawn)")
    else:
        p1 = outdir / "budget_curves.png"
        plot_curves(meta, rows, B, a, v, p1)
        print(f"wrote {p1}")

    p2 = outdir / "budget_bounds.png"
    plot_bounds(binary, a, v, p2)
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
