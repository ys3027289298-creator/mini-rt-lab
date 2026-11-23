#!/usr/bin/env python3
"""Convert an rt .raw frame to PNG (stdlib only). Usage: raw2png.py in.raw out.png [diffbase.raw]"""
import sys, struct, zlib

def read_raw(path):
    d = open(path, 'rb').read()
    if d[:4] != b'RT1\x00':
        raise SystemExit(f'{path}: not an RT1 frame')
    w, h, spp, hv = struct.unpack('<IIII', d[4:20])
    px = d[20:]
    if len(px) != w * h * 3:
        raise SystemExit(f'{path}: truncated ({len(px)} != {w*h*3})')
    return w, h, spp, hv, px

def write_png(path, w, h, rgb):
    raw = b''.join(b'\x00' + rgb[y*w*3:(y+1)*w*3] for y in range(h))
    def chunk(t, d):
        c = t + d
        return struct.pack('>I', len(d)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
    png = b'\x89PNG\r\n\x1a\n'
    png += chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 2, 0, 0, 0))
    png += chunk(b'IDAT', zlib.compress(raw, 9))
    png += chunk(b'IEND', b'')
    open(path, 'wb').write(png)

def main():
    if len(sys.argv) not in (3, 4):
        raise SystemExit(__doc__)
    w, h, spp, hv, px = read_raw(sys.argv[1])
    if len(sys.argv) == 4:
        w2, h2, _, _, base = read_raw(sys.argv[2])
        if (w, h) != (w2, h2):
            raise SystemExit('size mismatch')
        out = bytearray(len(px))
        for i in range(len(px)):
            d = abs(px[i] - base[i])
            out[i] = min(255, d * 4)
        px = bytes(out)
    write_png(sys.argv[-1], w, h, px)

if __name__ == '__main__':
    main()
