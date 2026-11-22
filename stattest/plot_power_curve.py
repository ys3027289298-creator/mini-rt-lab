#!/usr/bin/env python3
"""Render power_curve.csv (N delta power) to power_curve.png.

Line chart, 3 series (N=200/1000/5000), categorical slots 1-3 of the
reference palette (validated: ALL CHECKS PASS, light mode; aqua contrast
WARN relieved by direct labels on every line). Static PNG: the CSV is
the table view. DejaVu Sans stands in for the system sans.
"""
import csv
import sys

from PIL import Image, ImageDraw, ImageFont

W, H = 720, 440
ML, MR, MT, MB = 56, 96, 48, 48
PW, PH = W - ML - MR, H - MT - MB
X_MAX, Y_MAX = 0.30, 1.0

SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
SERIES = [("N=200", "#2a78d6"), ("N=1000", "#eb6834"), ("N=5000", "#1baf7a")]

FONT_DIR = "/usr/share/fonts/truetype/dejavu/"
F_TICK = ImageFont.truetype(FONT_DIR + "DejaVuSans.ttf", 11)
F_LABEL = ImageFont.truetype(FONT_DIR + "DejaVuSans.ttf", 11)
F_TITLE = ImageFont.truetype(FONT_DIR + "DejaVuSans-Bold.ttf", 14)


def sx(d):
    return ML + d / X_MAX * PW


def sy(p):
    return MT + (1 - p / Y_MAX) * PH


def main(path_in, path_out):
    data = {name: [] for name, _ in SERIES}
    with open(path_in) as f:
        for row in csv.reader(f, delimiter=" "):
            n, d, p = int(row[0]), float(row[1]), float(row[2])
            key = f"N={n}"
            if key in data:
                data[key].append((d, p))

    img = Image.new("RGB", (W, H), SURFACE)
    dr = ImageDraw.Draw(img)

    # title
    dr.text((ML, 14), "Power curve — chi-square goodness-of-fit (K=16, α=0.05)",
            font=F_TITLE, fill=INK)
    dr.text((ML, 32), "power vs chi-square divergence δ; dots mark δ_min at target power 0.80",
            font=F_TICK, fill=INK2)

    # gridlines + ticks
    for i in range(6):
        yv = i * 0.2
        y = sy(yv)
        dr.line([(ML, y), (ML + PW, y)], fill=GRID, width=1)
        dr.text((ML - 8, y), f"{yv:.1f}", font=F_TICK, fill=MUTED, anchor="rm")
    d = 0.0
    while d <= X_MAX + 1e-9:
        x = sx(d)
        dr.line([(x, MT), (x, MT + PH)], fill=GRID, width=1)
        dr.text((x, MT + PH + 8), f"{d:.2f}", font=F_TICK, fill=MUTED, anchor="ma")
        d += 0.05
    dr.line([(ML, MT + PH), (ML + PW, MT + PH)], fill=BASELINE, width=1)
    dr.text((ML + PW / 2, H - 10), "chi-square divergence δ",
            font=F_TICK, fill=MUTED, anchor="mb")
    dr.text((16, MT + PH / 2), "power", font=F_TICK, fill=MUTED, anchor="mm")

    # target power reference line
    yt = sy(0.8)
    dr.line([(ML, yt), (ML + PW, yt)], fill=BASELINE, width=1)
    dr.text((ML + PW - 4, yt - 5), "target power 0.80", font=F_TICK, fill=MUTED, anchor="rb")

    # series lines + end dots (2px surface ring)
    for name, color in SERIES:
        pts = [(sx(d), sy(p)) for d, p in data[name]]
        dr.line(pts, fill=color, width=2, joint="curve")
        xe, ye = pts[-1]
        dr.ellipse([xe - 6, ye - 6, xe + 6, ye + 6], fill=SURFACE)
        dr.ellipse([xe - 4, ye - 4, xe + 4, ye + 4], fill=color)

    # delta_min dots where each curve crosses power=0.8
    for name, color in SERIES:
        pts = data[name]
        for (d0, p0), (d1, p1) in zip(pts, pts[1:]):
            if p0 < 0.8 <= p1:
                dm = d0 + (0.8 - p0) * (d1 - d0) / (p1 - p0)
                x, y = sx(dm), yt
                dr.ellipse([x - 6, y - 6, x + 6, y + 6], fill=SURFACE)
                dr.ellipse([x - 4, y - 4, x + 4, y + 4], fill=color)
                break

    # direct labels at the x of max separation from other curves
    for name, color in SERIES:
        pts = data[name]
        others = [data[n] for n, _ in SERIES if n != name]
        i_best = max(range(len(pts)),
                     key=lambda i: min(abs(pts[i][1] - o[i][1]) for o in others))
        d, p = pts[i_best]
        below_others = all(p > o[i_best][1] for o in others)  # lower on screen
        lx, ly = sx(d), sy(p) + (14 if below_others else -14)
        dr.text((lx, ly), name, font=F_LABEL, fill=INK2,
                anchor="mb" if below_others else "ma")

    # legend (top-left inside plot)
    lx, ly = ML + 12, MT + 14
    for name, color in SERIES:
        dr.line([(lx, ly), (lx + 18, ly)], fill=color, width=2)
        dr.text((lx + 24, ly), name, font=F_LABEL, fill=INK2, anchor="lm")
        ly += 16

    img.save(path_out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "power_curve.csv",
         sys.argv[2] if len(sys.argv) > 2 else "power_curve.png")
