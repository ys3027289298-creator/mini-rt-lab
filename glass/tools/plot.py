#!/usr/bin/env python3
"""Render the moire experiment's CSVs into dependency-free SVG figures.

Usage: plot.py out/

Reads  out/summary.csv, out/stripe_D*_*.csv, out/samples_D*_*.csv
Writes out/stripe_D*.svg  — measured stripe profiles, grid vs offset, with the
                             predicted grid wavelength marked,
       out/samples_D*.svg — where the two samplers actually put their samples
                             relative to the panel edges, and
       out/summary.svg    — stripe spacing vs subdivision density, showing
                             lambda*D constant for the grid and scattered for
                             the offset.

Stdlib only; no numpy/matplotlib required.
"""
import csv
import os
import sys

GRID = "#c53030"
OFF = "#2b6cb0"


def read_profile(path):
    xs, ys = [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            xs.append(int(r["x"]))
            ys.append(float(r["residual"]))
    return xs, ys


def read_summary(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({k: float(v) for k, v in r.items()})
    return rows


def read_samples(path):
    """rows are (px, x, y) in absolute wall coordinates."""
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append((int(r["px"]), float(r["x"]), float(r["y"])))
    return rows


def samples_svg(panes, out_path, tier, pitch):
    """Show sample positions over a strip of wall, with panel edges drawn.

    panes: list of (label, rows). The point is visual: in the grid panel every
    pixel repeats one lattice, so samples hold a fixed offset from the panel
    edges and that offset drifts by frac(p) per pixel; in the offset panel it
    is shuffled pixel to pixel, breaking the lock between the two lattices.
    """
    margin_l, margin_t = 34, 66
    xmax = max(max(r[1] for r in rows) for _, rows in panes) + 1.0
    ymin = min(min(r[2] for r in rows) for _, rows in panes)
    ymax = max(max(r[2] for r in rows) for _, rows in panes) + 1.0
    scale = 20
    strip_h = (ymax - ymin) * scale
    gap = 60
    w = margin_l * 2 + len(panes) * (xmax * scale) + (len(panes) - 1) * gap
    h = margin_t + strip_h + 46
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}" font-family="sans-serif">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{w/2}" y="22" text-anchor="middle" font-size="16" '
        f'fill="#1a202c">sample positions on the wall — D={tier} '
        f'(panel pitch {pitch:.2f} px)</text>',
        f'<text x="{w/2}" y="42" text-anchor="middle" font-size="11" '
        f'fill="#666">dark vertical lines = panel edges, faint = pixel edges; '
        'dots = 16 samples per pixel, first 4 rows</text>',
    ]
    for pi, (label, rows) in enumerate(panes):
        pxl = margin_l + pi * (xmax * scale + gap)
        top = margin_t
        colour = GRID if pi == 0 else OFF
        parts.append(f'<text x="{pxl + xmax*scale/2:.0f}" y="{top - 10}" '
                     f'text-anchor="middle" font-size="13" fill="{colour}">'
                     f'{label}</text>')
        parts.append(f'<rect x="{pxl}" y="{top}" width="{xmax*scale}" '
                     f'height="{strip_h:.1f}" fill="#f7fafc" stroke="#cbd5e0"/>')
        edge = 0.0
        while edge <= xmax + 1e-9:
            ex = pxl + edge * scale
            parts.append(f'<line x1="{ex:.1f}" y1="{top}" x2="{ex:.1f}" '
                         f'y2="{top + strip_h:.1f}" stroke="#a0aec0"/>')
            edge += pitch
        for bx in range(0, int(xmax) + 1):
            ex = pxl + bx * scale
            parts.append(f'<line x1="{ex:.1f}" y1="{top}" x2="{ex:.1f}" '
                         f'y2="{top + strip_h:.1f}" stroke="#edf2f7"/>')
        for (px, x, y) in rows:
            cx = pxl + x * scale
            cy = top + (y - ymin) * scale
            parts.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="1.8" '
                         f'fill="{colour}" fill-opacity="0.7"/>')
    parts.append("</svg>")
    with open(out_path, "w") as f:
        f.write("\n".join(parts))


def profile_svg(tier, pred, band_grid, band_off, g, o, out_path):
    """Two stacked panels: the residual profile, and its amplitude bar."""
    w, h = 960, 380
    ml, mr, mt, mb = 60, 24, 58, 40
    pw = w - ml - mr
    ph = 190
    gx, gy = g
    ox, oy = o
    vals = gy + oy
    lo, hi = min(vals), max(vals)
    span = (hi - lo) or 1.0
    n = max(len(gx), len(ox))

    def X(i):
        return ml + pw * i / max(n - 1, 1)

    def Y(v, top):
        return top + ph - ph * (v - lo) / span

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}" font-family="sans-serif">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{ml}" y="22" font-size="15" fill="#1a202c">D={tier}: '
        f'column-averaged residual against the exact box filter</text>',
        f'<text x="{ml}" y="40" font-size="11" fill="#666">'
        f'predicted and measured grid stripe spacing = {pred:.2f} px</text>',
    ]
    top = mt
    parts.append(f'<text x="{ml}" y="{top - 8}" font-size="12" fill="{GRID}">'
                 f'grid (fixed lattice)</text>')
    pts = " ".join(f"{X(x):.1f},{Y(v, top):.1f}" for x, v in zip(gx, gy))
    parts.append(f'<polyline points="{pts}" fill="none" stroke="{GRID}" '
                 'stroke-width="1.2"/>')
    # mark the predicted wavelength with ticks along the top of the panel
    if pred > 0:
        x0 = ml
        xt = x0
        while xt < ml + pw:
            parts.append(f'<line x1="{xt:.1f}" y1="{top}" x2="{xt:.1f}" '
                         f'y2="{top + 12}" stroke="{GRID}" stroke-width="1" '
                         'stroke-dasharray="2 2"/>')
            xt += pred
        parts.append(f'<text x="{ml + pw}" y="{top + 10}" text-anchor="end" '
                     f'font-size="10" fill="{GRID}">ticks every {pred:.1f}px</text>')

    top2 = mt + ph + 64
    parts.append(f'<line x1="{ml}" y1="{mt + ph + 18}" x2="{w - mr}" '
                 f'y2="{mt + ph + 18}" stroke="#e2e8f0"/>')
    parts.append(f'<text x="{ml}" y="{top2 - 8}" font-size="12" fill="{OFF}">'
                 f'offset (R2 low-discrepancy), band RMS {band_off:.4f} vs '
                 f'{band_grid:.4f}</text>')
    pts = " ".join(f"{X(x):.1f},{Y(v, top2):.1f}" for x, v in zip(ox, oy))
    parts.append(f'<polyline points="{pts}" fill="none" stroke="{OFF}" '
                 'stroke-width="1.2"/>')
    parts.append(f'<text x="{ml}" y="{h - 12}" font-size="11" fill="#666">'
                 'x = wall column (px) →</text>')
    parts.append("</svg>")
    with open(out_path, "w") as f:
        f.write("\n".join(parts))


def summary_svg(rows, out_path):
    """lambda (y) vs D (x) on log-log, plus the lambda*D product column."""
    w, h = 880, 420
    ml, mr, mt, mb = 80, 170, 54, 60
    pw, ph = w - ml - mr, h - mt - mb
    Ds = [r["D"] for r in rows]
    lams = [r["grid_lam"] for r in rows] + [r["off_lam"] for r in rows]
    xlo, xhi = min(Ds) * 0.9, max(Ds) * 1.1
    ylo, yhi = min(lams) * 0.8, max(lams) * 1.25
    import math
    lx = lambda d: ml + pw * (math.log(d) - math.log(xlo)) / (math.log(xhi) - math.log(xlo))
    ly = lambda v: mt + ph - ph * (math.log(v) - math.log(ylo)) / (math.log(yhi) - math.log(ylo))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}" font-family="sans-serif">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{w/2}" y="24" text-anchor="middle" font-size="15" '
        'fill="#1a202c">stripe spacing vs subdivision density (log-log)</text>',
        f'<text x="{w/2}" y="42" text-anchor="middle" font-size="11" '
        'fill="#666">grid points lie on lambda ∝ 1/D; the offset has no '
        'coherent stripe, so its points are noise and mean nothing</text>',
    ]
    # axes
    parts.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#a0aec0"/>')
    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#a0aec0"/>')
    for r in rows:
        parts.append(f'<text x="{lx(r["D"]):.1f}" y="{mt+ph+18}" '
                     f'text-anchor="middle" font-size="11" fill="#4a5568">'
                     f'D={int(r["D"])}</text>')
    for frac in (0.0, 0.5, 1.0):
        v = math.exp(math.log(ylo) + frac * (math.log(yhi) - math.log(ylo)))
        y = ly(v)
        parts.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml+pw}" y2="{y:.1f}" stroke="#edf2f7"/>')
        parts.append(f'<text x="{ml-8}" y="{y+4:.1f}" text-anchor="end" '
                     f'font-size="11" fill="#4a5568">{v:.0f}</text>')
    # grid: fit lambda*D = C through the mean product and draw it
    gp = [r["grid_lamD"] for r in rows]
    cmean = sum(gp) / len(gp)
    pts = []
    xs = [xlo + (xhi - xlo) * i / 40 for i in range(41)]
    for d in xs:
        lam = cmean / d
        if ylo <= lam <= yhi:
            pts.append(f"{lx(d):.1f},{ly(lam):.1f}")
    parts.append(f'<polyline points="{" ".join(pts)}" fill="none" '
                 f'stroke="{GRID}" stroke-width="1.2" stroke-dasharray="5 4"/>')
    for r in rows:
        parts.append(f'<circle cx="{lx(r["D"]):.1f}" cy="{ly(r["grid_lam"]):.1f}" '
                     f'r="5" fill="{GRID}"/>')
        parts.append(f'<circle cx="{lx(r["D"]):.1f}" cy="{ly(r["off_lam"]):.1f}" '
                     f'r="5" fill="{OFF}"/>')
    parts.append(f'<text x="{ml+8}" y="{mt+14}" font-size="12" fill="{GRID}">'
                 f'grid (λ·D ≈ {cmean:.0f}, spread '
                 f'{max(gp)/min(gp):.2f}×)</text>')
    op = [r["off_lamD"] for r in rows]
    parts.append(f'<text x="{ml+8}" y="{mt+32}" font-size="12" fill="{OFF}">'
                 f'offset (no stripe; λ·D scatter {max(op)/min(op):.2f}×)</text>')
    bo = [r["band_off"] for r in rows]
    bg = [r["band_grid"] for r in rows]
    parts.append(f'<text x="{ml+8}" y="{mt+50}" font-size="12" fill="#2d3748">'
                 f'band RMS: grid {min(bg):.3f}–{max(bg):.3f}, offset '
                 f'{min(bo):.3f}–{max(bo):.3f}</text>')
    # side table of lambda*D
    tx = ml + pw + 20
    parts.append(f'<text x="{tx}" y="{mt+8}" font-size="11" fill="#2d3748">λ·D</text>')
    for i, r in enumerate(rows):
        yy = mt + 26 + i * 20
        parts.append(f'<text x="{tx}" y="{yy}" font-size="10" fill="{GRID}">'
                     f'{r["grid_lamD"]:.0f}</text>')
        parts.append(f'<text x="{tx+54}" y="{yy}" font-size="10" fill="{OFF}">'
                     f'{r["off_lamD"]:.0f}</text>')
    parts.append("</svg>")
    with open(out_path, "w") as f:
        f.write("\n".join(parts))


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "out"
    summary = os.path.join(out, "summary.csv")
    if not os.path.exists(summary):
        print("missing", summary)
        return 1
    rows = read_summary(summary)
    made = 0
    for r in rows:
        tier = int(r["D"])
        pg = os.path.join(out, f"stripe_D{tier}_grid.csv")
        po = os.path.join(out, f"stripe_D{tier}_offset.csv")
        if os.path.exists(pg) and os.path.exists(po):
            profile_svg(tier, r["pred_lam"], r["band_grid"], r["band_off"],
                        read_profile(pg), read_profile(po),
                        os.path.join(out, f"stripe_D{tier}.svg"))
            made += 1
        sg = os.path.join(out, f"samples_D{tier}_grid.csv")
        so = os.path.join(out, f"samples_D{tier}_offset.csv")
        if os.path.exists(sg) and os.path.exists(so):
            samples_svg([("grid: fixed lattice", read_samples(sg)),
                         ("offset: R2 low-discrepancy", read_samples(so))],
                        os.path.join(out, f"samples_D{tier}.svg"), tier,
                        r["pitch"])
    summary_svg(rows, os.path.join(out, "summary.svg"))
    print(f"wrote {made} stripe profiles + summary.svg to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
