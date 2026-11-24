#!/usr/bin/env python3
"""Edge-block statistics for the nested-sphere renders.

The walk-mode render with the rho-derived limits is the reference answer for
the same scene (identical geometry, camera and per-pixel RNG layout). Every
other render is compared against it inside the reference silhouette.

Two block metrics, both anchored on the reference rather than on hand-picked
luminance thresholds:

  gap pixels      reference-silhouette pixels the test render never reached
                  (its hit mask is 0 there). These are the hard chunks cut out
                  of the inner sphere's edge; the legacy one-query tracer
                  loses whole pixels this way.
  diff-block px   pixels whose luminance differs from the reference by more
                  than 3 * sigma_ref, where sigma_ref is the Monte-Carlo noise
                  floor measured inside the reference interior (with a 2
                  code-value floor so a noiseless render cannot false-fire).
                  This catches partial energy loss that gap counting misses
                  at high SPP, where a pixel is reached by only a few samples.

Usage: analyze.py OUT_DIR REF.ppm TEST1.ppm [TEST2.ppm ...]
Reads the matching .mask.ppm next to every color PPM. Writes edge crops with
gap pixels marked red / diff-block pixels marked orange, plus stats.json and
stats.txt.
"""
import json
import os
import subprocess
import sys
from statistics import median


def read_ppm(path):
    with open(path, "rb") as f:
        if f.readline().strip() != b"P6":
            raise ValueError(f"{path}: not a P6 PPM")
        line = f.readline()
        while line.startswith(b"#"):
            line = f.readline()
        w, h = map(int, line.split())
        if int(f.readline()) != 255:
            raise ValueError(f"{path}: unsupported maxval")
        buf = f.read()
    if len(buf) != w * h * 3:
        raise ValueError(f"{path}: size mismatch")
    return w, h, buf


def lum(buf, i):
    r, g, b = buf[3 * i], buf[3 * i + 1], buf[3 * i + 2]
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def neighbors(i, w, h):
    x, y = i % w, i // w
    if x > 0:
        yield i - 1
    if x + 1 < w:
        yield i + 1
    if y > 0:
        yield i - w
    if y + 1 < h:
        yield i + w


def erode(mask, w, h, iters):
    m = bytearray(mask)
    for _ in range(iters):
        out = bytearray(len(m))
        for i in range(len(m)):
            if m[i] and all(m[j] for j in neighbors(i, w, h)):
                out[i] = 1
        m = out
    return m


def write_ppm(path, w, h, buf):
    with open(path, "wb") as f:
        f.write(f"P6\n{w} {h}\n255\n".encode())
        f.write(bytes(buf))


def to_png(ppm_path):
    png_path = ppm_path[:-4] + ".png"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i",
                    ppm_path, png_path], check=True)


def main(out_dir, ref_path, test_paths):
    os.makedirs(out_dir, exist_ok=True)
    w, h, ref_buf = read_ppm(ref_path)
    _, _, ref_mbuf = read_ppm(ref_path[:-4] + ".mask.ppm")
    n = w * h
    ref_mask = bytearray(1 if ref_mbuf[3 * i] > 127 else 0 for i in range(n))
    L_ref = [lum(ref_buf, i) for i in range(n)]

    interior = erode(ref_mask, w, h, 2)
    edge_band = bytearray(ref_mask[i] and not interior[i] for i in range(n))

    # Reference MC noise floor from neighbor differences in its own interior.
    diffs = []
    for i in range(n):
        if interior[i]:
            for j in neighbors(i, w, h):
                if j > i and interior[j]:
                    diffs.append(L_ref[i] - L_ref[j])
    mean_d = sum(diffs) / len(diffs)
    var = sum((d - mean_d) ** 2 for d in diffs) / len(diffs)
    sigma_ref = (var ** 0.5) / (2 ** 0.5)
    diff_thr = max(3.0 * sigma_ref, 2.0)
    L_inner = median(L_ref[i] for i in range(n) if interior[i])

    # Crop window around the reference silhouette.
    xs = [i % w for i in range(n) if ref_mask[i]]
    ys = [i // w for i in range(n) if ref_mask[i]]
    margin = 14
    x0, x1 = max(0, min(xs) - margin), min(w, max(xs) + margin + 1)
    y0, y1 = max(0, min(ys) - margin), min(h, max(ys) + margin + 1)
    cw, ch = x1 - x0, y1 - y0

    def crop_image(buf, gap_set, diff_set, path):
        crop = bytearray(3 * cw * ch)
        for yy in range(ch):
            for xx in range(cw):
                si = (yy + y0) * w + (xx + x0)
                di = 3 * (yy * cw + xx)
                crop[di:di + 3] = buf[3 * si:3 * si + 3]
        for i in diff_set:                    # orange: luminance-level blocks
            xx, yy = i % w - x0, i // w - y0
            if 0 <= xx < cw and 0 <= yy < ch:
                di = 3 * (yy * cw + xx)
                crop[di], crop[di + 1], crop[di + 2] = 255, 140, 0
        for i in gap_set:                     # red: fully missing pixels
            xx, yy = i % w - x0, i // w - y0
            if 0 <= xx < cw and 0 <= yy < ch:
                di = 3 * (yy * cw + xx)
                crop[di], crop[di + 1], crop[di + 2] = 255, 0, 0
        write_ppm(path, cw, ch, crop)
        to_png(path)
        os.remove(path)

    stats = []
    jobs = [(os.path.basename(ref_path)[:-4], ref_buf, ref_mask, True)]
    for p in test_paths:
        tw, th, tbuf = read_ppm(p)
        assert (tw, th) == (w, h), f"{p}: resolution mismatch"
        _, _, tmb = read_ppm(p[:-4] + ".mask.ppm")
        tmask = bytearray(1 if tmb[3 * i] > 127 else 0 for i in range(n))
        jobs.append((os.path.basename(p)[:-4], tbuf, tmask, False))

    edge_idx = [i for i in range(n) if edge_band[i]]
    for name, buf, mask, is_ref in jobs:
        gap = {i for i in range(n) if ref_mask[i] and not mask[i]}
        gap_edge = {i for i in gap if edge_band[i]}
        if is_ref:
            diff_blocks = set()
            edge_mean_d = edge_rms_d = sil_mean_d = sil_rms_d = 0.0
        else:
            L = [lum(buf, i) for i in range(n)]
            diff_blocks = {i for i in range(n)
                           if ref_mask[i] and abs(L[i] - L_ref[i]) > diff_thr}
            # Bias metrics: systematic energy loss shifts every pixel the same
            # way, so the 3-sigma block count under-reports it; mean/RMS of
            # the luminance difference measure the loss directly.
            ed = [L[i] - L_ref[i] for i in edge_idx]
            edge_mean_d = sum(ed) / len(ed)
            edge_rms_d = (sum(v * v for v in ed) / len(ed)) ** 0.5
            sd = [L[i] - L_ref[i] for i in range(n) if ref_mask[i]]
            sil_mean_d = sum(sd) / len(sd)
            sil_rms_d = (sum(v * v for v in sd) / len(sd)) ** 0.5
        diff_edge = {i for i in diff_blocks if edge_band[i]}
        crop_image(buf, gap_edge, diff_edge,
                   os.path.join(out_dir, name + ".edge.ppm"))
        stats.append({"render": name,
                      "is_reference": is_ref,
                      "mask_px": sum(mask),
                      "gap_px": len(gap),
                      "gap_px_edge": len(gap_edge),
                      "diff_block_px": len(diff_blocks),
                      "diff_block_px_edge": len(diff_edge),
                      "edge_mean_dL": round(edge_mean_d, 4),
                      "edge_rms_dL": round(edge_rms_d, 4),
                      "sil_mean_dL": round(sil_mean_d, 4),
                      "sil_rms_dL": round(sil_rms_d, 4)})

    meta = {"crop": [x0, y0, x1, y1],
            "sigma_ref": round(sigma_ref, 4),
            "diff_threshold": round(diff_thr, 3),
            "L_inner_ref": round(L_inner, 3),
            "renders": stats}
    with open(os.path.join(out_dir, "stats.json"), "w") as f:
        json.dump(meta, f, indent=2)

    hdr = ("render".ljust(30) + "gap@edge".rjust(9) + "dblk@edge".rjust(10) +
           "mean dL".rjust(9) + "rms dL".rjust(9) +
           "sil mean".rjust(9) + "sil rms".rjust(9))
    lines = [f"sigma_ref={sigma_ref:.4f}  diff_threshold(3sigma,min2)={diff_thr:.3f}",
             "(dL in 8-bit code values, measured on the edge band / full silhouette)",
             hdr, "-" * len(hdr)]
    for s in stats:
        lines.append(s["render"].ljust(30) +
                     str(s["gap_px_edge"]).rjust(9) +
                     str(s["diff_block_px_edge"]).rjust(10) +
                     f"{s['edge_mean_dL']:.3f}".rjust(9) +
                     f"{s['edge_rms_dL']:.3f}".rjust(9) +
                     f"{s['sil_mean_dL']:.3f}".rjust(9) +
                     f"{s['sil_rms_dL']:.3f}".rjust(9))
    table = "\n".join(lines)
    print(table)
    with open(os.path.join(out_dir, "stats.txt"), "w") as f:
        f.write(table + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit("usage: analyze.py OUT_DIR REF.ppm TEST.ppm [TEST.ppm ...]")
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
