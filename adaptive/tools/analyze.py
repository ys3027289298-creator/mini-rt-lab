#!/usr/bin/env python3
"""Summarise adaptive_pixels.csv: early-stop fractions, stop_n distribution,
bias/diff quantiles and tolerance breaches.  Standard library only; run with
    python3 -I tools/analyze.py [artifacts/adaptive_pixels.csv]
"""
import csv
import sys
from collections import defaultdict


def pct(xs, q):
    """Empirical percentile q in [0,100], linear interpolation."""
    if not xs:
        return float("nan")
    xs = sorted(xs)
    pos = (len(xs) - 1) * q / 100.0
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(xs):
        return xs[-1]
    return xs[lo] * (1 - frac) + xs[lo + 1] * frac


def main(path):
    groups = defaultdict(list)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            groups[(row["scene"], row["criterion"])].append(row)

    hdr = ("scene", "crit", "n", "early%", "stop_n p25/p50/p75/max",
           "bias p50/p95/max", "diff p95/max", "viol_truth%", "viol_full%")
    print("{:<6} {:<6} {:>5} {:>7} {:>22} {:>20} {:>16} {:>11} {:>11}".format(*hdr))
    for (scene, crit), rows in sorted(groups.items()):
        n = len(rows)
        early = sum(int(r["early"]) for r in rows)
        stop_n = [int(r["stop_n"]) for r in rows]
        bias = [float(r["bias"]) for r in rows]
        diff = [float(r["diff"]) for r in rows]
        vt = sum(int(r["viol_truth"]) for r in rows)
        vf = sum(int(r["viol_full"]) for r in rows)
        print("{:<6} {:<6} {:>5} {:>6.1f}% {:>22} {:>20} {:>16} {:>10.1f}% {:>10.1f}%".format(
            scene, crit, n, 100.0 * early / n,
            "{:.0f}/{:.0f}/{:.0f}/{:.0f}".format(
                pct(stop_n, 25), pct(stop_n, 50), pct(stop_n, 75), max(stop_n)),
            "{:.4f}/{:.4f}/{:.4f}".format(pct(bias, 50), pct(bias, 95), max(bias)),
            "{:.4f}/{:.4f}".format(pct(diff, 95), max(diff)),
            100.0 * vt / n, 100.0 * vf / n))

    # paired EB-vs-naive comparison on the same pixels
    print("\npaired per-pixel comparison (naive relative to eb):")
    by_pixel = defaultdict(dict)
    for (scene, crit), rows in groups.items():
        for r in rows:
            by_pixel[(scene, int(r["pixel"]))][crit] = r
    for scene in sorted({k[0] for k in by_pixel}):
        worse = saved = total = 0
        for (sc, _), pair in by_pixel.items():
            if sc != scene or "eb" not in pair or "naive" not in pair:
                continue
            total += 1
            if float(pair["naive"]["bias"]) > float(pair["eb"]["bias"]):
                worse += 1
            saved += int(pair["eb"]["stop_n"]) - int(pair["naive"]["stop_n"])
        print(f"  {scene}: naive worse on {worse}/{total} pixels; "
              f"eb uses {saved} extra samples total vs naive")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "artifacts/adaptive_pixels.csv")
