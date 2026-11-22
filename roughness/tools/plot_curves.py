#!/usr/bin/env python3
"""Plot the perceptual-roughness -> GGX-alpha mapping and its validation.

Reads the CSVs emitted by the C self-test (run `make test` first) and writes
out/curves.png and out/validation.png. Palette follows the project dataviz
convention: blue/orange/aqua categorical slots, single-hue sequential fills.
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# the C program writes out/ relative to the working directory, which is the
# module root when run via `make test`; fall back to a root-relative path.
OUT = os.path.join(HERE, "..", "out")
if not os.path.isdir(OUT):
    OUT = os.path.join(os.getcwd(), "out")

# validated categorical slots (light mode): blue, orange, aqua
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["system-ui", "Segoe UI", "DejaVu Sans", "Arial"],
    "font.size": 11,
    "text.color": INK,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def load(name):
    with open(os.path.join(OUT, name), newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def main():
    sweep = load("alpha_sweep.csv")
    curve = load("mapping_curve.csv")
    val = load("validation.csv")
    cross = load("render_crosscheck.csv")
    summary = open(os.path.join(OUT, "summary.txt")).read()
    p_line = next(l for l in summary.splitlines() if "r^" in l)
    p_str = p_line.split("r^")[1].split()[0].rstrip("),")

    # ---- figure 1: mapping curves + measured blur response --------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2))

    ax1.plot(curve["r"], curve["alpha_lut"], color=BLUE, lw=2.2,
             label="LUT mapping (shipped)")
    ax1.plot(curve["r"], curve["alpha_power"], color=ORANGE, lw=2,
             ls="--", label=f"power approx  $\\alpha=r^{{{p_str}}}$")
    ax1.plot(curve["r"], curve["alpha_r2"], color=AQUA, lw=2, ls=":",
             label="textbook  $\\alpha=r^2$")
    ax1.plot([0, 1], [0, 1], color=MUTED, lw=1, alpha=0.6)
    ax1.annotate("r = 0 → perfect mirror", xy=(0, 0), xytext=(0.08, 0.42),
                 color=MUTED, fontsize=9,
                 arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
    ax1.set_title("Perceptual roughness $r$ → GGX $\\alpha$")
    ax1.set_xlabel("perceptual roughness  r")
    ax1.set_ylabel("GGX alpha")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.legend(frameon=False, fontsize=9.5, loc="upper left")

    ax2.plot(sweep["alpha"], sweep["blur_normalized"], color=BLUE, lw=2.2)
    ax2.set_title("Measured blur response  B($\\alpha$)")
    ax2.set_xlabel("GGX alpha")
    ax2.set_ylabel("normalized NDF median angle")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1.02)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "curves.png"), dpi=150)
    plt.close(fig)

    # ---- figure 2: validation errors + rendered cross-check -------------
    fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5.2))

    ax3.axhline(0, color=MUTED, lw=1)
    ax3.plot(val["r"], [100 * e for e in val["error_lut"]], color=BLUE,
             lw=2.2, marker="o", ms=5, label="LUT mapping")
    ax3.plot(val["r"], [100 * e for e in val["error_power"]], color=ORANGE,
             lw=2, ls="--", marker="s", ms=5, label=f"power $r^{{{p_str}}}$")
    ax3.plot(val["r"], [100 * e for e in val["error_r2"]], color=AQUA,
             lw=2, ls=":", marker="^", ms=5, label="textbook $r^2$")
    ax3.set_title("Measured blur error vs target (lower is better)")
    ax3.set_xlabel("perceptual roughness  r")
    ax3.set_ylabel("measured blur − target  (%)")
    ax3.legend(frameon=False, fontsize=9.5)

    ax4.plot(cross["alpha"], cross["pixel_sigma_px"], color=ORANGE, lw=2.2,
             marker="o", ms=5, label="rendered pixel σ")
    ax4.set_title("MC-rendered highlight σ (cross-check)")
    ax4.set_xlabel("GGX alpha")
    ax4.set_ylabel("pixel-space σ (px)")
    ax4.annotate("peak from Smith G2 masking\n(physical, not the metric)",
                 xy=(0.5, max(cross["pixel_sigma_px"])),
                 xytext=(0.12, max(cross["pixel_sigma_px"]) * 0.78),
                 color=MUTED, fontsize=9,
                 arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
    ax4.legend(frameon=False, fontsize=9.5)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "validation.png"), dpi=150)
    plt.close(fig)

    print("wrote", os.path.join(OUT, "curves.png"))
    print("wrote", os.path.join(OUT, "validation.png"))


if __name__ == "__main__":
    sys.exit(main())
