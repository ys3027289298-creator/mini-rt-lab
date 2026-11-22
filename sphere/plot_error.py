#!/usr/bin/env python3
"""Plot sphere-intersection error data produced by src/sphere.c.

Usage: python3 plot_error.py [csv_dir] [out_dir]
Defaults: ../artifacts/sphere for both.

Figures:
  fig1_sweep.png   absolute root error vs ray magnitude, two geometries
  fig2_tangent.png near-tangent behaviour around D = 0, three cancellation scales
  fig3_random.png  binned relative-error envelopes, stable vs naive, 2e5 rays
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# validated categorical palette (dataviz skill, all-pairs gate, 3-slot cap)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE, GRID, INK, MUTED, GRAY = "#fcfcfb", "#e1e0d9", "#0b0b0b", "#52514e", "#898781"
EPS = np.float32(np.finfo(np.float32).eps)

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "grid.color": GRID, "grid.linewidth": 0.8,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Liberation Sans", "Arial"],
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
})


def load(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None,
                         encoding="utf-8", loose=True)


# ---------------- Fig 1: magnitude sweeps ----------------
def fig_sweep(d, out):
    sw = load(os.path.join(d, "sweep.csv"))
    # Above s ~ 2^12 the float ray direction can no longer encode the grazing
    # offset (eps*s*s overwhelms the O(1) geometry), and the double reference
    # alternates hit/miss on grid rounding. Only the reference-reliable range
    # is plotted; rows beyond it carry status 3 in the CSV.
    S_CUT = 4096.0
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))

    for ax, geom, title, slopes in [
        (axes[0], "grazing", "grazing ray, h = 0.5 R  (true D = O(1), b² and c ~ s²)",
         [(2.0, INK, "slope ε·s²  (discriminant rounding)"),
          (1.0, AQUA, "slope √ε·s  (tangent-band clamp)")]),
        (axes[1], "central", "central ray, h = 0  (true D ≈ 1, well conditioned)",
         [(1.0, INK, "slope ε·s  (ordinary rounding)")]),
    ]:
        m = (sw["geom"] == geom) & (sw["s"].astype(float) <= S_CUT)
        s = sw["s"][m].astype(float)
        stable = sw["abs_stable"][m].astype(float)
        naive = sw["abs_naive"][m].astype(float)
        bound = sw["bound"][m].astype(float)
        status = sw["status"][m].astype(int)
        two = status == 1
        tan = status == 2
        nv = ~np.isnan(naive)
        nlost = np.isnan(naive)

        ax.loglog(s[two], stable[two], "o", ms=4.5, color=BLUE,
                  label="stable: two distinct roots")
        ax.loglog(s[tan], stable[tan], "s", ms=5, color=AQUA,
                  label="stable: merged tangent (band)")
        ax.loglog(s[nv], naive[nv], "^", ms=4.5, color=ORANGE,
                  label="naive −b ± √D")
        if nlost.any():
            # mark lost hits on the axis floor, not at a fabricated y value
            ax.loglog(s[nlost],
                      np.full(nlost.sum(), np.nanmin(stable) * 0.25),
                      "x", ms=10, mew=2.2, color=ORANGE,
                      label="naive: hit lost (D computed < 0)")
        ax.loglog(s, bound, "--", lw=1.6, color=GRAY,
                  label="analytic bound  Kε(scale/√D + |t| + S), K = 64")

        sgrid = np.logspace(np.log10(s.min()), np.log10(s.max()), 40)
        for slope, col, lab in slopes:
            if slope == 2.0:
                ref = EPS * sgrid**2
            elif slope == 1.0 and geom == "grazing":
                ref = np.sqrt(EPS) * sgrid          # √ε·s clamp scale
            else:
                ref = EPS * sgrid
            ax.loglog(sgrid, ref, ":", lw=1.3, color=col, label=lab)

        ax.set_xlim(s.min() * 0.7, S_CUT * 1.4)
        ax.text(S_CUT * 1.04, 3e-11,
                "s > 4096 omitted:\nfloat ray direction\nitself loses the\ngrazing angle",
                ha="left", va="bottom", color=MUTED, fontsize=8.3)

        ax.set_xlabel("ray magnitude s   (origin distance, R = 1)")
        ax.set_ylabel("absolute error of near root  |t − t_ref|")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.6)
        ax.legend(fontsize=7.6, framealpha=0.94, loc="upper left")
    fig.suptitle("Root error vs t magnitude: cancellation drives the grazing "
                 "error as ε·s² until the tangent band caps it at √ε·s",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1_sweep.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------- Fig 2: near-tangent band ----------------
def fig_tangent(d, out):
    tg = load(os.path.join(d, "tangent.csv"))
    Bs = [1.0, 1e3, 1e6]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.0))
    for ax, B in zip(axes, Bs):
        m = tg["B"] == B
        Drel = tg["Drel"][m].astype(float)
        side = tg["side"][m].astype(int)
        status = tg["status"][m].astype(int)
        abs_near = tg["abs_near"][m].astype(float)
        poly = tg["poly_res"][m].astype(float)
        scale = B * B + 1.0
        tol = 8 * EPS * scale

        pos_two = (side == 0) & (status == 1)
        pos_tan = (side == 0) & (status == 2)
        neg_tan = (side == 1) & (status == 2)
        neg_miss = (side == 1) & (status == 0)

        ax.loglog(Drel[pos_two], abs_near[pos_two], "o", ms=5, color=BLUE,
                  label="+ side: two roots, |Δt| vs reference")
        ax.loglog(Drel[pos_tan], abs_near[pos_tan], "s", ms=5.5, color=AQUA,
                  label="+ side: merged tangent")
        floor = 1e-12
        ax.loglog(Drel[neg_tan], np.maximum(poly[neg_tan], floor), "x",
                  ms=8, mew=2, color=ORANGE,
                  label="− side clamp: |F(t)| at reported root")
        if neg_miss.any():
            ax.loglog(Drel[neg_miss], np.full(neg_miss.sum(), floor), "|",
                      ms=8, mew=1.8, color=GRAY, label="− side: clear miss")

        Dg = np.logspace(-9, 0, 80)
        ax.loglog(Dg, 64 * EPS * scale / np.sqrt(Dg), "--", lw=1.5, color=GRAY,
                  label="bound  Kε·scale/√D")
        ax.axhline(np.sqrt(EPS * scale), ls=":", lw=1.4, color=ORANGE,
                   label="√(ε·scale): clamp displacement scale")
        ax.axvspan(Drel.min(), min(tol, 1.0), color=ORANGE, alpha=0.07)
        if tol <= 1.5:
            ax.axvline(tol, color=ORANGE, lw=1.2, alpha=0.9)
            ax.text(tol * 1.3, 2e-11, f"band edge  8ε·scale = {tol:.1e}",
                    color="#c24e1d", fontsize=8.2, va="bottom")
        else:
            ax.text(0.98, 0.30,
                    f"band edge 8ε·scale = {tol:.1e} > 1:\nentire sampled range "
                    "is numerically\ntangent — no root is resolvable",
                    transform=ax.transAxes, ha="right", va="top",
                    color="#c24e1d", fontsize=8.2,
                    bbox=dict(facecolor=SURFACE, edgecolor="#e8c9bb", boxstyle="round,pad=0.35"))

        ax.set_xlim(Drel.min() * 0.7, 2)
        ax.set_ylim(5e-13, 2e3)
        ax.set_xlabel("|D|   (magnitude of the discriminant)")
        ax.set_title(f"B = {B:g}\nscale ≈ B² = {scale:.2e}", pad=8)
        ax.grid(True, which="both", alpha=0.6)
        ax.legend(fontsize=7.0, loc="upper right", framealpha=0.94)
    axes[0].set_ylabel("root error  |Δt|   /   polynomial residual  |F(t)|")
    fig.suptitle("Near-tangent policy: |D| ≤ 8ε·scale merges to one tangent root "
                 "on BOTH sides of D = 0; beyond the band, D < 0 stays a miss",
                 fontsize=12, y=1.04)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(out, "fig2_tangent.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------- Fig 3: randomized rays ----------------
def fig_random(d, out):
    rd = load(os.path.join(d, "random.csv"))
    t = np.abs(rd["t_ref"].astype(float))
    rs = rd["rel_stable"].astype(float)
    rn = rd["rel_naive"].astype(float)

    edges = np.logspace(np.log10(t.min()), np.log10(t.max()), 9)
    idx = np.digitize(t, edges)
    centers = np.sqrt(edges[:-1] * edges[1:])

    def envelope(vals, q):
        out = np.full(len(centers), np.nan)
        for i in range(1, len(edges)):
            sel = vals[idx == i]
            sel = sel[np.isfinite(sel)]
            if sel.size:
                out[i - 1] = np.percentile(sel, q)
        return out

    fig, ax = plt.subplots(figsize=(11, 5.2))
    ax.loglog(centers, envelope(rs, 99), "o-", color=BLUE, ms=5, lw=1.6,
              label="stable solver, 99th percentile")
    ax.loglog(centers, envelope(rs, 100), "s--", color=BLUE, ms=5, lw=1.2,
              alpha=0.65, label="stable solver, worst case")
    ax.loglog(centers, envelope(rn, 99), "o-", color=ORANGE, ms=5, lw=1.6,
              label="naive −b ± √D, 99th percentile")
    ax.loglog(centers, envelope(rn, 100), "^--", color=ORANGE, ms=5, lw=1.2,
              alpha=0.65, label="naive −b ± √D, worst case")
    ax.axhline(EPS, color=INK, ls=":", lw=1.3,
               label="ε ≈ 1.19e-7: relative floor reached at large |t|")
    ax.axvspan(t.min(), 0.05, color=GRAY, alpha=0.10)
    ax.text(0.012, 4e-4,
            "small-|t| regime:\nnormalisation floor ε·S dominates,\nso the "
            "guarantee is the ABSOLUTE bound\nof Fig. 1 (Kε·scale/√D), not a "
            "relative ε",
            ha="left", va="top", color=MUTED, fontsize=8.6)
    ax.set_xlabel("|t_ref|   (near-root hit distance)")
    ax.set_ylabel("normalised error  |t − t_ref| / max(|t_ref|, ε·S)")
    ax.set_title("181,020 reference-hit random rays, magnitudes 2⁻⁴…2²⁰: on generic "
                 "rays the two solvers are statistically tied (median ratio 1.0); "
                 "stable earns its cost on grazing/tangent geometry, where naive "
                 "loses the 1 hit seen in Fig. 1",
                 fontsize=11.5)
    ax.grid(True, which="both", alpha=0.6)
    ax.legend(fontsize=8.2, framealpha=0.94, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3_random.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    csv_dir = sys.argv[1] if len(sys.argv) > 1 else "../artifacts/sphere"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else csv_dir
    fig_sweep(csv_dir, out_dir)
    fig_tangent(csv_dir, out_dir)
    fig_random(csv_dir, out_dir)
    print("wrote fig1_sweep.png fig2_tangent.png fig3_random.png to", out_dir)
