#!/usr/bin/env python3
"""Render the diff histogram CSV into a dependency-free SVG.

Usage: histogram.py diff_hist.csv out/diff_hist.svg
Stdlib only; no numpy/matplotlib required.
"""
import csv
import sys


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 1
    bins, fracs, labels = [], [], []
    with open(sys.argv[1], newline="") as f:
        for row in csv.DictReader(f):
            bins.append(int(row["bin"]))
            fracs.append(float(row["fraction"]))
            labels.append(f"{float(row['lo']):.0f}-{float(row['hi']):.0f}")

    w, h = 900, 460
    ml, mr, mt, mb = 70, 20, 40, 90
    pw, ph = w - ml - mr, h - mt - mb
    peak = max(fracs) or 1.0
    bw = pw / len(bins)
    cumsum = 0.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}" font-family="sans-serif">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{w/2}" y="24" text-anchor="middle" font-size="18" '
        'fill="#222">per-pixel |delta luminance| histogram (sRGB levels, '
        "path-trace vs rasterize)</text>",
    ]
    # y grid + percent labels
    for g in range(5):
        yv = peak * g / 4
        y = mt + ph - ph * g / 4
        parts.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{w-mr}" y2="{y:.1f}" '
                     'stroke="#ddd"/>')
        parts.append(f'<text x="{ml-8}" y="{y+4:.1f}" text-anchor="end" '
                     f'font-size="11" fill="#555">{yv*100:.1f}%</text>')
    # bars, colored by cumulative mass: first 90% of pixels green-ish, tail red
    for i, (b, fr) in enumerate(zip(bins, fracs)):
        cumsum += fr
        bh = ph * fr / peak
        x = ml + i * bw
        y = mt + ph - bh
        color = "#4c9f70" if cumsum <= 0.90 else ("#e0a458" if cumsum <= 0.99 else "#c0504d")
        parts.append(f'<rect x="{x+1:.1f}" y="{y:.1f}" width="{bw-2:.1f}" '
                     f'height="{bh:.1f}" fill="{color}"/>')
        if i % 2 == 0:
            parts.append(f'<text x="{x+bw/2:.1f}" y="{h-mb+16}" '
                         f'text-anchor="middle" font-size="9" fill="#555" '
                         f'transform="rotate(60 {x+bw/2:.1f},{h-mb+16})">'
                         f"{labels[i]}</text>")
    # cumulative line
    cumsum = 0.0
    pts = []
    for i, fr in enumerate(fracs):
        cumsum += fr
        x = ml + (i + 0.5) * bw
        y = mt + ph - ph * cumsum
        pts.append(f"{x:.1f},{y:.1f}")
    parts.append(f'<polyline points="{" ".join(pts)}" fill="none" '
                 'stroke="#2b6cb0" stroke-width="1.6"/>')
    parts.append(f'<text x="{w-mr}" y="{mt+12}" text-anchor="end" font-size="11" '
                 'fill="#2b6cb0">blue line: cumulative fraction</text>')
    parts.append(f'<text x="{ml+ph*0}" y="{h-18}" font-size="12" fill="#333">'
                 "green = first 90% of pixel mass, amber = 90-99%, red = tail"
                 "</text>")
    parts.append("</svg>")
    with open(sys.argv[2], "w") as f:
        f.write("\n".join(parts))
    print(f"wrote {sys.argv[2]} ({len(bins)} bins, peak {peak*100:.2f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
