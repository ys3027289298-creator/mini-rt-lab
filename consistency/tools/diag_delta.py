#!/usr/bin/env python3
"""Diagnostic: distribution of per-pixel color/depth deltas between frames.

Usage: python3 tools/diag_delta.py rast.fbuf trace.fbuf
"""

import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from consist_check import read_fbuf


def main():
    r = read_fbuf(sys.argv[1])
    t = read_fbuf(sys.argv[2])
    rid = r[:, :, 4].round().astype(int)
    tid = t[:, :, 4].round().astype(int)
    same = (rid == tid) & (rid >= 0)
    d = np.abs(r[:, :, 0:3].astype(float) - t[:, :, 0:3].astype(float)).max(axis=2)
    ds = d[same]
    print("same-id covered pixels:", int(same.sum()))
    if ds.size:
        for q in (0.0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0):
            print(f"q={q:5.2f}  color_delta={np.quantile(ds, q):.6f}")
    for tri in np.unique(rid[same]):
        m = same & (rid == tri)
        dd = d[m]
        print(f"tri {tri}: n={int(m.sum()):6d} min={dd.min():.6f} "
              f"max={dd.max():.6f} std={dd.std():.6f}")
    z = np.abs(r[:, :, 3].astype(float) - t[:, :, 3].astype(float))[same]
    if z.size:
        print("depth abs delta max/mean:", z.max(), z.mean())


if __name__ == "__main__":
    main()
