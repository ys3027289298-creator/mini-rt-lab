"""
Figures for the shadow-map re-allocation experiment.

Reads shadowmap/out/*.csv and writes shadowmap/figs/*.png.

Palette: dataviz reference instance, slots 1-3 (validated all-pairs, light mode).
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

OUT = os.path.join(os.path.dirname(__file__), "..", "out")
FIG = os.path.join(os.path.dirname(__file__), "..", "figs")
os.makedirs(FIG, exist_ok=True)

# --- palette (light mode) -------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
S1 = "#2a78d6"   # blue
S2 = "#eb6834"   # orange
S3 = "#1baf7a"   # aqua
S4 = "#eda100"   # yellow (single-series context only)
S7 = "#6250d6"   # violet

rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 10,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def read_csv(name):
    with open(os.path.join(OUT, name), newline="") as f:
        return list(csv.DictReader(f))


def style(ax, title, sub=None, xlabel=None, ylabel=None, logx=False, logy=False):
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="major", alpha=0.9)
    ax.grid(True, which="minor", alpha=0.35, linewidth=0.5)
    nlines = sub.count("\n") + 1 if sub else 0
    ax.set_title(title, color=INK, fontsize=12, loc="left",
                 pad=12 + 13 * nlines)
    if sub:
        ax.text(0, 1.012, sub, transform=ax.transAxes, color=INK2, fontsize=9,
                va="bottom", ha="left", linespacing=1.35)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.tick_params(which="both", length=3)
    ax.set_axisbelow(True)


def header(fig, title, sub):
    """Title + wrapped subtitle, positioned from real text heights (inches)."""
    hin = fig.get_figheight()
    t_in = 13 * 1.35 / 72.0                 # title block
    n = sub.count("\n") + 1
    s_in = 9 * 1.35 * n / 72.0              # subtitle block
    y_title = 0.975
    fig.text(0.012, y_title, title, ha="left", va="top",
             color=INK, fontsize=13)
    # subtitle top sits 1.15 title-lines below the title top
    y_sub = y_title - (13 * 1.15 + 2) / 72.0 / hin
    fig.text(0.012, y_sub, sub, ha="left", va="top",
             color=INK2, fontsize=9, linespacing=1.4)
    return y_sub - s_in / hin - 0.03


# ---------------------------------------------------------------- figure 1
def fig_exponent():
    """Derived local exponent alpha(r) against the constants tested."""
    rows = read_csv("density.csv")
    r = [float(x["r"]) for x in rows]
    al = [float(x["alpha_local"]) for x in rows]
    fit = float(rows[0]["alpha_bandfit"])

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.plot(r, al, color=S1, lw=2, label="derived  α(r) = 3r²/(r²+H²)")
    ax.axhline(fit, color=S2, lw=2, ls=(0, (5, 2)),
               label=f"constant band fit  α = {fit:.2f}")
    for a, lab in [(1.25, "α = 1.25"), (2.0, "α = 2.0")]:
        ax.axhline(a, color=MUTED, lw=1, ls=(0, (2, 3)))
        ax.text(2.46, a, lab, color=MUTED, fontsize=8, va="center", ha="right",
                bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1))
    ax.axvline(1.0, color=AXIS, lw=0.8)
    ax.text(1.03, 2.75, "r = H", color=MUTED, fontsize=8)
    ax.annotate("α → 0 at the light", xy=(0.03, 0.02),
                xytext=(0.035, 2.2), color=INK2, fontsize=8,
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8))
    style(ax, "The exponent is not constant",
          "Local exponent of the equal-solid-angle allocation, H = 1.\n"
          "Only a sweep can decide which constant is best.",
          "wall radius r / H", "texel-area exponent  α", logx=True)
    ax.set_ylim(-0.15, 3.15)
    ax.set_xlim(0.009, 2.6)
    ax.legend(loc="upper left", labelcolor=INK2, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_exponent.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- figure 2
def fig_density():
    """Derived, band-fit and observed texel density on one axis."""
    d = read_csv("density.csv")
    obs = read_csv("density_obs.csv")
    r = [float(x["r"]) for x in d]
    nd = [float(x["n_derived"]) for x in d]
    nb = [float(x["n_bandfit"]) for x in d]

    by = {}
    for x in obs:
        by.setdefault(x["tag"], []).append((float(x["r"]), float(x["coverage"])))

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.plot(r, nd, color=S1, lw=2,
            label="derived density  n ∝ r / (r²+H²)$^{3/2}$")
    ax.plot(r, nb, color=S2, lw=2, ls=(0, (5, 2)),
            label="constant band fit  α = 0.36")
    ax.plot(*zip(*by["derived_exact"]), color=S1, lw=0, marker="o", ms=5,
            mfc=SURFACE, mew=1.4, label="measured, exact derived rows")
    ax.plot(*zip(*by["planar_a0"]), color=S3, lw=0, marker="o", ms=5,
            mfc=SURFACE, mew=1.4, label="measured, uniform-area rows")

    ax.axvline(0.16, color=AXIS, lw=0.8)
    ax.text(0.175, 5.2e4, "rows thinner\nthan one texel",
            color=MUTED, fontsize=8, va="top")
    style(ax, "Predicted and measured texel density agree",
          "Texels per unit wall area. Measured points are the row allocation read\n"
          "at the 20 sample radii; the scatter below r ≈ 0.16H is one-row quantisation.",
          "wall radius r / H", "texels per wall area  n  (H$^{-2}$)", logx=True)
    ax.set_ylim(1.4e3, 6.2e4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
              labelcolor=INK2, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_density.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- figure 3
def fig_exponent_sweep():
    """Artifact area vs allocation exponent, for the three intensities."""
    rows = read_csv("sweep.csv")
    a = [float(x["alpha"]) for x in rows]
    near = [float(x["leak_near"]) for x in rows]
    far = [float(x["leak_far"]) for x in rows]
    tot = [float(x["leak_H2"]) for x in rows]
    fit = float(read_csv("density.csv")[0]["alpha_bandfit"])

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.3), sharey=False)
    specs = [
        (near, "near the light  (r < 0.5H)", S1,
         "the patch the report is about"),
        (far, "far from the light  (r ≥ 0.5H)", S2,
         "shadows out here lose texels"),
        (tot, "whole wall", S7, "the two terms together"),
    ]
    for ax, (y, title, col, note) in zip(axes, specs):
        ax.plot(a, y, color=col, lw=2)
        ax.plot(a, y, color=col, lw=0, marker="o", ms=4,
                mfc=SURFACE, mew=1.2)
        i = min(range(len(y)), key=lambda k: y[k])
        ax.plot([a[i]], [y[i]], marker="o", ms=9, color=col, mfc=col,
                mec=SURFACE, mew=2, zorder=5)
        ax.annotate(f"α* = {a[i]:.2f}\n{y[i]:.4f}",
                    xy=(a[i], y[i]), xytext=(8, 12),
                    textcoords="offset points", color=INK2, fontsize=8,
                    arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8))
        ax.axvline(fit, color=MUTED, lw=1, ls=(0, (2, 3)))
        style(ax, title, None, "allocation exponent  α",
              "leak area  (H²)" if ax is axes[0] else None)
        ax.set_ylim(0, max(y) * 1.12)
    axes[0].text(0.03, 0.30, "α=0 uniform · α>0 near-dense · α<0 reversed",
                 transform=axes[0].transAxes, color=MUTED, fontsize=7.5,
                 va="top")
    top = header(fig,
                 "Where the minimum sits depends on which region is counted",
                 "Leak = wall area a lookup leaves lit although the analytic shadow is dark,\n"
                 "averaged over 12 randomised blocker layouts. Dotted line: the band fit α = 0.36.")
    fig.tight_layout(rect=(0, 0, 1, top))
    fig.savefig(os.path.join(FIG, "fig3_sweep.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- figure 4
def fig_intensity():
    """Visible bright-patch area at three intensities, named allocations."""
    rows = read_csv("points.csv")
    tags = ["planar_a0", "derived_exact", "const_bandfit", "const_a1.25",
            "const_a2.0", "reverse_am1"]
    nice = {
        "planar_a0": "uniform α=0",
        "derived_exact": "derived (exact)",
        "const_bandfit": "constant α=0.36",
        "const_a1.25": "constant α=1.25",
        "const_a2.0": "constant α=2.0",
        "reverse_am1": "reversed α=−1",
    }
    ints = ["low", "mid", "high"]
    data = {t: {i: [] for i in ints} for t in tags}
    for x in rows:
        data[x["tag"]][x["intensity"]].append(
            (float(x["r_mid"]), float(x["vis_dens"])))

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.2), sharey=True)
    for ax, inten, sub in zip(axes, ints,
                              ["I = 0.35", "I = 1.0", "I = 3.0"]):
        style_map = {
            "derived_exact": dict(color=S1, lw=2),
            "planar_a0": dict(color=S2, lw=2),
            "const_bandfit": dict(color=S3, lw=1.4, ls=(0, (5, 2))),
            "const_a1.25": dict(color=S3, lw=1.4, ls=(0, (1, 2))),
            "const_a2.0": dict(color=S3, lw=1.4, ls=(0, (3, 1, 1, 1))),
        }
        for t in tags:
            if t == "reverse_am1":
                continue
            xs, ys = zip(*data[t][inten])
            ax.plot(xs, ys, **style_map[t])
        xs, ys = zip(*data["reverse_am1"][inten])
        ax.plot(xs, ys, color=S7, lw=1.4, ls=(0, (4, 2)))
        style(ax, sub, None, "wall radius r / H",
              "visible leak per wall area  (H⁻¹)" if ax is axes[0] else None,
              logx=True)
    # one legend for the whole figure, identity via line + label
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], color=S1, lw=2, label="derived, exact warp"),
        Line2D([], [], color=S2, lw=2, label="uniform α=0"),
        Line2D([], [], color=S3, lw=1.4, ls=(0, (5, 2)), label="constant α=0.36"),
        Line2D([], [], color=S3, lw=1.4, ls=(0, (1, 2)), label="constant α=1.25"),
        Line2D([], [], color=S3, lw=1.4, ls=(0, (3, 1, 1, 1)), label="constant α=2.0"),
        Line2D([], [], color=S7, lw=1.4, ls=(0, (4, 2)), label="reversed α=−1"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=6,
               bbox_to_anchor=(0.5, 0.0), labelcolor=INK2, fontsize=8)
    top = header(fig, "Near the light the derived warp removes the bright ring; "
                      "far out it pays for it",
                 "Visible leak density = per-ring leak (lit where truth is shadowed, and dimmed\n"
                 "radiance above 0.10) divided by the ring area, so near and far rings compare.")
    fig.tight_layout(rect=(0, 0.07, 1, top))
    fig.savefig(os.path.join(FIG, "fig4_intensity.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- figure 6
def fig_pareto():
    """Near cost against far cost, one point per exponent."""
    rows = read_csv("sweep.csv")
    a = [float(x["alpha"]) for x in rows]
    near = [float(x["leak_near"]) / 0.25 for x in rows]
    far = [float(x["leak_far"]) / 1.75 for x in rows]
    fit = float(read_csv("density.csv")[0]["alpha_bandfit"])

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#0d366b"])
    sc = ax.scatter(near, far, c=a, cmap=cmap, s=52, edgecolors=SURFACE,
                    linewidths=1.5, zorder=4)
    ax.plot(near, far, color=GRID, lw=1, zorder=1)
    for i in (0, 3, 6, 8, 10, 12, 14):
        ax.annotate(f"α={a[i]:.2f}", (near[i], far[i]),
                    textcoords="offset points",
                    xytext=(7, 10) if i != 6 else (7, -13),
                    color=INK2, fontsize=8)
    # circle the exponent that minimises each side separately
    ax.set_ylim(top=3.2e-1)
    for vals, col, lab, off in ((near, S1, "best near", (28, -52)),
                                (far, S2, "best far", (56, -80))):
        j = min(range(len(vals)), key=lambda k: vals[k])
        ax.plot([near[j]], [far[j]], marker="o", ms=14, mfc="none",
                mec=col, mew=2, zorder=6)
        ax.annotate(f"{lab}  α={a[j]:.2f}", (near[j], far[j]),
                    textcoords="offset points", xytext=off,
                    color=col, fontsize=8.5, ha="left",
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.8,
                                    connectionstyle="angle3,angleA=0,angleB=90"))
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("allocation exponent α", color=INK2)
    style(ax, "No single exponent is best on both sides",
          "Both axes are leak density (leak area ÷ region area), log scale.\n"
          "Lower-left is better. The two circles sit at different α.",
          "near-light leak density  (r < 0.5H)", "far-field leak density  (r ≥ 0.5H)",
          logx=True, logy=True)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig6_pareto.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- figure 5
def fig_summary_bars():
    """Plain comparison of the six allocations, near band."""
    rows = read_csv("summary.csv")
    tags = ["planar_a0", "derived_exact", "const_bandfit", "const_a1.25",
            "const_a2.0", "reverse_am1"]
    nice = {
        "planar_a0": "uniform α=0",
        "derived_exact": "derived, exact",
        "const_bandfit": "constant α=0.36",
        "const_a1.25": "constant α=1.25",
        "const_a2.0": "constant α=2.0",
        "reverse_am1": "reversed α=−1",
    }
    by = {x["tag"]: x for x in rows}
    near = [float(by[t]["leak_near"]) for t in tags]
    far = [float(by[t]["leak_far"]) for t in tags]
    sd = [float(by[t]["leak_near_sd"]) for t in tags]
    lab = [nice[t] for t in tags]

    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    x = range(len(tags))
    w = 0.36
    ax.bar([i - w / 2 for i in x], far, width=w, color=S2, label="far  (r ≥ 0.5H)")
    ax.bar([i + w / 2 for i in x], near, width=w, color=S1,
           label="near the light  (r < 0.5H)")
    ax.errorbar([i + w / 2 for i in x], near, yerr=sd, fmt="none",
                ecolor=INK2, elinewidth=1, capsize=3)
    for i, v in enumerate(near):
        ax.text(i + w / 2, v + sd[i] + 0.0006, f"{v:.4f}", ha="center",
                va="bottom", fontsize=7.5, color=INK2)
    ax.set_xticks(list(x))
    ax.set_xticklabels(lab, fontsize=9)
    style(ax, "The near-light patch and the far-field cost trade off",
          "Leak area per allocation. Error bars: s.d. over 12 randomised layouts.\n"
          "α = 1.25 is the best compromise; the exact warp is best near the light.",
          None, "leak area  (H²)")
    ax.legend(loc="upper left", labelcolor=INK2, fontsize=9)
    ax.set_ylim(0, max(far) * 1.18)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_summary.png"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    fig_exponent()
    fig_density()
    fig_exponent_sweep()
    fig_intensity()
    fig_pareto()
    fig_summary_bars()
    print("figures:", sorted(os.listdir(FIG)))
