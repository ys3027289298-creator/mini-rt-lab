#!/usr/bin/env python3
"""render_check.py — 纯标准库位图渲染（zlib 写 PNG），用于无浏览器环境的目视检查。

功能：
  1. 把 pipeline 产出的 PPM 中间图转成 PNG；
  2. 按 analyze_pipeline 的同一数据，把逐级亮度曲线栅格化成 PNG
     （自带 5x7 点阵字体，不依赖任何字体文件或第三方库）。

用法：
  python3 -I render_check.py <run_dir> [<run_dir> ...]
输出：每个 run 目录下 cut*.png，以及 runs 目录旁的 brightness_curve.png
"""
import os
import struct
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_pipeline import (  # noqa: E402
    analyze_run, stage_deltas, read_pfm, STAGES, DROP_THRESHOLD,
)

# ---------------- 极简 PNG 编码 ----------------

def write_png(path, w, h, rows):
    """rows: 长度 h 的列表，每项是长度 w*3 的 bytes（RGB）。"""
    def chunk(typ, data):
        c = typ + data
        return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c))
    raw = b"".join(b"\x00" + bytes(r) for r in rows)
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(chunk(b"IHDR", ihdr))
        f.write(chunk(b"IDAT", zlib.compress(raw, 9)))
        f.write(chunk(b"IEND", b""))


def read_ppm(path):
    with open(path, "rb") as f:
        assert f.readline().strip() == b"P6"
        line = f.readline()
        while line.startswith(b"#"):
            line = f.readline()
        w, h = map(int, line.split())
        assert int(f.readline()) == 255
        data = f.read(w * h * 3)
    return w, h, data


def ppm_to_png(ppm_path, png_path, scale=2):
    w, h, data = read_ppm(ppm_path)
    rows = []
    for y in range(h):
        row = bytearray()
        for x in range(w):
            i = (y * w + x) * 3
            row += data[i:i + 3] * scale
        band = bytes(row)
        for _ in range(scale):
            rows.append(band)
    write_png(png_path, w * scale, h * scale, rows)


# ---------------- 5x7 点阵字体（仅本图需要的字符） ----------------

FONT = {
    "0": ["01110", "10001", "10011", "10101", "11001", "10001", "01110"],
    "1": ["00100", "01100", "00100", "00100", "00100", "00100", "01110"],
    "2": ["01110", "10001", "00001", "00010", "00100", "01000", "11111"],
    "3": ["11110", "00001", "00001", "01110", "00001", "00001", "11110"],
    "4": ["00010", "00110", "01010", "10010", "11111", "00010", "00010"],
    "5": ["11111", "10000", "11110", "00001", "00001", "10001", "01110"],
    "6": ["00110", "01000", "10000", "11110", "10001", "10001", "01110"],
    "7": ["11111", "00001", "00010", "00100", "01000", "01000", "01000"],
    "8": ["01110", "10001", "10001", "01110", "10001", "10001", "01110"],
    "9": ["01110", "10001", "10001", "01111", "00001", "00010", "01100"],
    ".": ["00000", "00000", "00000", "00000", "00000", "00110", "00110"],
    "-": ["00000", "00000", "00000", "11111", "00000", "00000", "00000"],
    "%": ["11001", "11010", "00100", "00100", "01000", "10110", "10011"],
    "c": ["00000", "00000", "01110", "10001", "10000", "10001", "01110"],
    "u": ["00000", "00000", "10001", "10001", "10001", "10011", "01101"],
    "t": ["01000", "01000", "11100", "01000", "01000", "01010", "00100"],
    " ": ["00000"] * 7,
    "a": ["00000", "00000", "01110", "00001", "01111", "10001", "01111"],
    "b": ["10000", "10000", "10110", "11001", "10001", "10001", "11110"],
    "c": ["00000", "00000", "01110", "10001", "10000", "10001", "01110"],
    "d": ["00001", "00001", "01101", "10011", "10001", "10011", "01101"],
    "e": ["00000", "00000", "01110", "10001", "11111", "10000", "01110"],
    "f": ["00110", "01001", "01000", "11100", "01000", "01000", "01000"],
    "g": ["00000", "01111", "10001", "10001", "01111", "00001", "01110"],
    "h": ["10000", "10000", "10110", "11001", "10001", "10001", "10001"],
    "i": ["00100", "00000", "01100", "00100", "00100", "00100", "01110"],
    "l": ["01100", "00100", "00100", "00100", "00100", "00100", "01110"],
    "m": ["00000", "00000", "11010", "10101", "10101", "10001", "10001"],
    "n": ["00000", "00000", "10110", "11001", "10001", "10001", "10001"],
    "o": ["00000", "00000", "01110", "10001", "10001", "10001", "01110"],
    "p": ["00000", "11110", "10001", "10001", "11110", "10000", "10000"],
    "r": ["00000", "00000", "10111", "11001", "10000", "10000", "10000"],
    "s": ["00000", "00000", "01111", "10000", "01110", "00001", "11110"],
    "v": ["00000", "00000", "10001", "10001", "10001", "01010", "00100"],
    "x": ["00000", "00000", "10001", "01010", "00100", "01010", "10001"],
    "y": ["00000", "10001", "10001", "01111", "00001", "00010", "01100"],
}


def draw_text(buf, w, h, x, y, text, color, scale=1):
    for ci, ch in enumerate(text):
        glyph = FONT.get(ch, FONT[" "])
        for gy, row in enumerate(glyph):
            for gx, bit in enumerate(row):
                if bit == "1":
                    for sy in range(scale):
                        for sx in range(scale):
                            px, py = x + ci * 6 * scale + gx * scale + sx, y + gy * scale + sy
                            if 0 <= px < w and 0 <= py < h:
                                buf[py][px * 3:px * 3 + 3] = bytes(color)


def draw_line(buf, w, h, x0, y0, x1, y1, color):
    n = max(abs(x1 - x0), abs(y1 - y0), 1)
    for i in range(n + 1):
        x, y = round(x0 + (x1 - x0) * i / n), round(y0 + (y1 - y0) * i / n)
        if 0 <= x < w and 0 <= y < h:
            buf[y][x * 3:x * 3 + 3] = bytes(color)


def draw_circle(buf, w, h, cx, cy, r, color, ring=None):
    for y in range(cy - r - 1, cy + r + 2):
        for x in range(cx - r - 1, cx + r + 2):
            if 0 <= x < w and 0 <= y < h and (x - cx) ** 2 + (y - cy) ** 2 <= r * r:
                c = ring if ring and (x - cx) ** 2 + (y - cy) ** 2 > (r - 2) ** 2 else color
                buf[y][x * 3:x * 3 + 3] = bytes(c)


def hex_rgb(s):
    return [int(s[i:i + 2], 16) for i in (1, 3, 5)]


# ---------------- 曲线栅格化 ----------------

def render_curve(runs, out_path):
    W, H = 900, 460
    ml, mr, mt, mb = 92, 30, 40, 78
    pw, ph = W - ml - mr, H - mt - mb
    ink, ink2, grid, axis = (11, 11, 11), (82, 81, 78), (225, 224, 217), (195, 194, 183)
    surface = (252, 252, 251)
    danger = hex_rgb("#e34948")
    colors = [hex_rgb(c) for c in ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")]

    buf = [bytearray(surface) * W for _ in range(H)]
    all_y = [p[2] for r in runs for p in r["points"]]
    ymax = max(all_y) * 1.12

    def x_of(i):
        return ml + round(pw * i / 4.0)

    def y_of(v):
        return mt + round(ph * (1 - v / ymax))

    for g in range(5):  # 网格与 y 刻度
        v = ymax * g / 4.0
        yy = y_of(v)
        draw_line(buf, W, H, ml, yy, ml + pw, yy, grid)
        draw_text(buf, W, H, ml - 62, yy - 3, f"{v:.4f}", ink2)
    draw_line(buf, W, H, ml, mt, ml, mt + ph, axis)
    draw_line(buf, W, H, ml, mt + ph, ml + pw, mt + ph, axis)
    for i in range(5):  # x 轴 cut 标号
        xx = x_of(i)
        draw_text(buf, W, H, xx - 10, mt + ph + 14, f"cut{i+1}", ink, scale=2)
        draw_text(buf, W, H, xx - 4 * len(STAGES[i]) - 2, mt + ph + 38, STAGES[i], ink2)

    for ri, run in enumerate(runs):
        color = colors[ri % len(colors)]
        coords = [(x_of(p[0]), y_of(p[2])) for p in run["points"]]
        for a, b in zip(coords, coords[1:]):
            draw_line(buf, W, H, a[0], a[1], b[0], b[1], color)
        for (x, y), p in zip(coords, run["points"]):
            draw_circle(buf, W, H, x, y, 5, color, ring=surface)
        # 数值直接标注（墨色文字放数据点下方，不靠颜色辨认系列；
        # 多系列重合时按系列号垂直错开）
        for (x, y), p in zip(coords, run["points"]):
            draw_text(buf, W, H, x - 22, y + 10 + ri * 11, f"{p[2]:.4f}", ink)
        lx, ly = coords[-1]
        draw_text(buf, W, H, lx - 8 * len(run["name"]) - 4, ly + 26 + ri * 11, run["name"], ink)

    for row in stage_deltas(runs[0]["points"]):  # 超阈值落差标记
        stage_i, prev, cur, drop, bad = row
        if bad:
            yb = y_of(cur)
            draw_text(buf, W, H, x_of(stage_i) - 14, yb - 26, f"-{drop*100:.1f}%", danger, scale=2)

    write_png(out_path, W, H, [bytes(r) for r in buf])


def main():
    run_dirs = sys.argv[1:]
    runs = []
    for d in run_dirs:
        runs.append({"dir": d, "name": os.path.basename(os.path.normpath(d)),
                     "points": analyze_run(d)})
        for cut in range(1, 6):
            ppm = os.path.join(d, f"cut{cut}.ppm")
            png = os.path.join(d, f"cut{cut}.png")
            ppm_to_png(ppm, png, scale=3)
            print(f"{ppm} -> {png}")
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.normpath(run_dirs[0]))),
                       "brightness_curve.png")
    render_curve(runs, out)
    print(f"curve -> {out}")


if __name__ == "__main__":
    main()
