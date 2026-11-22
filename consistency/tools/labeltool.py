#!/usr/bin/env python3
"""Labeling aid: dump candidate diff pixels per class with diagnostics.

It never writes the label file itself. The workflow is:
  1. render both implementations at the label's base resolution
  2. run this tool, read the candidate coordinates and diagnostics
  3. hand-verify each candidate (visually / against the scene geometry)
  4. write the curated labels/*.json file with provenance filled in

Keeping the human in the loop is what makes labels/*.json a ground-truth
set rather than a copy of the automatic classifier's output.
"""

import argparse
import json

import numpy as np

from consist_check import classify, read_fbuf, CLASSES


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("rast_fbuf")
    ap.add_argument("trace_fbuf")
    ap.add_argument("--per-class", type=int, default=12,
                    help="max candidates printed per class")
    args = ap.parse_args()

    rast = read_fbuf(args.rast_fbuf)
    trace = read_fbuf(args.trace_fbuf)
    h, w = rast.shape[0], rast.shape[1]
    labels, metrics = classify(rast, trace)

    print(f"# resolution {w}x{h}")
    for cls in CLASSES:
        ys, xs = np.nonzero(labels == cls)
        print(f"# {cls}: {len(ys)} pixels")
        if len(ys) == 0:
            continue
        # spread samples across the image instead of taking a clump
        order = np.linspace(0, len(ys) - 1, min(args.per_class, len(ys))).round().astype(int)
        for i in order:
            y, x = int(ys[i]), int(xs[i])
            print(json.dumps({
                "u": round(x / (w - 1), 6),
                "v": round(y / (h - 1), 6),
                "x": x, "y": y,
                "cls": cls,
                "id_rast": int(rast[y, x, 4]),
                "id_trace": int(trace[y, x, 4]),
                "color_delta": round(float(metrics["color_delta"][y, x]), 6),
                "depth_rel_delta": round(float(metrics["depth_delta"][y, x]), 6),
                "on_boundary": bool(metrics["boundary"][y, x]),
            }, ensure_ascii=False))


if __name__ == "__main__":
    main()
