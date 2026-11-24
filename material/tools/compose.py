#!/usr/bin/env python3
"""Compose side-by-side highlight crops with variance labels.

Reads highlight_old.pfm / highlight_new.pfm + stats.json, writes:
  crop_old.png, crop_new.png   tone-mapped highlight crops (2x)
  highlight_compare.svg        both crops side by side, variance labels

Usage: compose.py <artifacts_dir>
"""
import base64
import json
import struct
import sys
import zlib


def read_pfm(path):
    with open(path, "rb") as f:
        assert f.readline().strip() == b"PF"
        w, h = map(int, f.readline().split())
        scale = float(f.readline())
        assert scale < 0  # little-endian
        data = f.read()
    n = w * h * 3
    vals = struct.unpack("<%df" % n, data[: 4 * n])
    # PFM rows are bottom-up; keep only the R channel (grayscale scene)
    img = [[vals[((h - 1 - y) * w + x) * 3] for x in range(w)] for y in range(h)]
    return w, h, img


def tonemap(v):
    v = max(v, 0.0)
    v = v / (1.0 + v)
    return int(round(255 * v ** (1 / 2.2)))


def write_png(path, w, h, rgb):
    raw = b"".join(b"\x00" + bytes(rgb[y * w * 3:(y + 1) * w * 3]) for y in range(h))

    def chunk(tag, data):
        c = tag + data
        return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c))

    png = (b"\x89PNG\r\n\x1a\n"
           + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
           + chunk(b"IDAT", zlib.compress(raw, 9))
           + chunk(b"IEND", b""))
    with open(path, "wb") as f:
        f.write(png)
    return png


def crop_scaled(img, cx, cy, cw, ch, scale):
    out = []
    for y in range(ch * scale):
        for x in range(cw * scale):
            g = tonemap(img[cy + y // scale][cx + x // scale])
            out.extend((g, g, g))
    return out


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "../artifacts"
    with open(f"{outdir}/stats.json") as f:
        stats = json.load(f)
    c = stats["crop"]
    w, h, img_old = read_pfm(f"{outdir}/highlight_old.pfm")
    _, _, img_new = read_pfm(f"{outdir}/highlight_new.pfm")

    scale = 2
    cw, ch = c["w"] * scale, c["h"] * scale
    rgb_old = crop_scaled(img_old, c["x"], c["y"], c["w"], c["h"], scale)
    rgb_new = crop_scaled(img_new, c["x"], c["y"], c["w"], c["h"], scale)
    png_old = write_png(f"{outdir}/crop_old.png", cw, ch, rgb_old)
    png_new = write_png(f"{outdir}/crop_new.png", cw, ch, rgb_new)

    g01 = next(g for g in stats["groups"] if g["rough"] == 0.1)
    v_old, v_new = g01["var_cos"], g01["var_vndf_tb"]
    ratio = v_old / v_new if v_new > 0 else float("inf")

    b64_old = base64.b64encode(png_old).decode()
    b64_new = base64.b64encode(png_new).decode()
    pad, label_h, gap = 12, 46, 24
    Wtot = 2 * cw + 3 * pad + gap
    Htot = ch + 2 * pad + label_h
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
     width="{Wtot}" height="{Htot}" viewBox="0 0 {Wtot} {Htot}">
  <rect width="100%" height="100%" fill="#1a1a1a"/>
  <image x="{pad}" y="{pad}" width="{cw}" height="{ch}"
         xlink:href="data:image/png;base64,{b64_old}"/>
  <image x="{2 * pad + cw + gap}" y="{pad}" width="{cw}" height="{ch}"
         xlink:href="data:image/png;base64,{b64_new}"/>
  <text x="{pad}" y="{ch + pad + 22}" fill="#eee" font-family="monospace" font-size="14">old: cosine sampling</text>
  <text x="{pad}" y="{ch + pad + 40}" fill="#fc0" font-family="monospace" font-size="14">noise var = {v_old:.3e}</text>
  <text x="{2 * pad + cw + gap}" y="{ch + pad + 22}" fill="#eee" font-family="monospace" font-size="14">new: VNDF + Taylor grazing band</text>
  <text x="{2 * pad + cw + gap}" y="{ch + pad + 40}" fill="#fc0" font-family="monospace" font-size="14">noise var = {v_new:.3e}  ({ratio:.1f}x lower)</text>
</svg>
"""
    with open(f"{outdir}/highlight_compare.svg", "w") as f:
        f.write(svg)
    print(f"wrote {outdir}/crop_old.png, crop_new.png, highlight_compare.svg")
    print(f"crop variance: old={v_old:.4e} new={v_new:.4e} ratio={ratio:.2f}")


if __name__ == "__main__":
    main()
