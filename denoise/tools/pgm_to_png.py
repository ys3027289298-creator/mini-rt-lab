#!/usr/bin/env python3
"""PGM(P5) -> PNG，仅用标准库，方便直接查看噪点图像。

用法: python3 -I tools/pgm_to_png.py in.pgm out.png
"""
import struct
import sys
import zlib


def read_pgm(path):
    with open(path, "rb") as f:
        assert f.readline().strip() == b"P5"
        line = f.readline()
        while line.startswith(b"#"):
            line = f.readline()
        w, h = map(int, line.split())
        maxv = int(f.readline())
        data = f.read()
    assert maxv == 255 and len(data) == w * h
    return w, h, data


def write_png(path, w, h, gray):
    def chunk(tag, payload):
        c = tag + payload
        return struct.pack(">I", len(payload)) + c + struct.pack(">I", zlib.crc32(c) & 0xffffffff)

    raw = bytearray()
    for y in range(h):  # 每行前加 filter type 0
        raw.append(0)
        raw.extend(gray[y * w:(y + 1) * w])
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0)  # 8-bit grayscale
    png = (b"\x89PNG\r\n\x1a\n"
           + chunk(b"IHDR", ihdr)
           + chunk(b"IDAT", zlib.compress(bytes(raw), 9))
           + chunk(b"IEND", b""))
    with open(path, "wb") as f:
        f.write(png)


if __name__ == "__main__":
    w, h, data = read_pgm(sys.argv[1])
    write_png(sys.argv[2], w, h, data)
    print(f"wrote {sys.argv[2]} {w}x{h}")
