#!/usr/bin/env python3
"""
Cross-check the three acceptance artifacts against the shared scene.

Inputs (all must exist under <artdir>, or the run fails):
  histogram.csv     sampling module  -> histogram
  variance.json     material module  -> variance
  image.ppm         geometry module  -> rendered image
  image_hash.json   geometry module  -> FNV-1a hash of image.ppm

Checks:
  1. all four files present (missing one => failure, per completion rule)
  2. every artifact references the same scene_id and seed; stream ids are
     1/2/3 (sampling/material/geometry), proving one scene, three streams
  3. histogram: bin counts sum to the declared draw count, uniform spread
  4. variance: within tolerance of the analytic Var(cos)=1/18
  5. image hash: recompute FNV-1a64 over the exact PPM bytes and compare

Runs under `python3 -I`; stdlib only.
"""
import json
import os
import sys


def fnv1a64(path):
    h = 14695981039346656037
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            for byte in chunk:
                h ^= byte
                h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h


def fail(msg, failures):
    failures.append(msg)
    print(f"FAIL: {msg}")


def main(artdir):
    failures = []
    expected_files = ["histogram.csv", "variance.json",
                      "image.ppm", "image_hash.json"]
    paths = {n: os.path.join(artdir, n) for n in expected_files}

    print("== presence ==")
    for n, p in paths.items():
        ok = os.path.isfile(p)
        print(f"{n}: {'present' if ok else 'MISSING'}")
        if not ok:
            fail(f"missing artifact {n}", failures)
    if failures:
        # Nothing else is trustworthy if a dataset is absent.
        print("\nRESULT: FAIL (a required dataset is missing)")
        return 1

    print("\n== shared scene linkage ==")
    # histogram.csv carries key=value header comments
    hist_meta = {}
    with open(paths["histogram.csv"]) as f:
        for line in f:
            if not line.startswith("#"):
                break
            for kv in line.lstrip("# ").strip().split(","):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    hist_meta[k.strip()] = v.strip()
    var_meta = json.load(open(paths["variance.json"]))
    img_meta = json.load(open(paths["image_hash.json"]))

    scene_ids = {
        "histogram": hist_meta.get("scene_id"),
        "variance": var_meta.get("scene_id"),
        "image_hash": img_meta.get("scene_id"),
    }
    print("scene_ids:", scene_ids)
    if len(set(scene_ids.values())) != 1:
        fail("scene_id differs across artifacts", failures)

    seeds = {
        "histogram": int(hist_meta.get("seed", "-1")),
        "variance": int(var_meta.get("seed", -1)),
        "image_hash": int(img_meta.get("seed", -1)),
    }
    print("seeds:", seeds)
    if len(set(seeds.values())) != 1:
        fail("seed differs across artifacts", failures)

    streams = {
        "histogram": int(hist_meta.get("stream", -1)),
        "variance": int(var_meta.get("stream", -1)),
        "image_hash": int(img_meta.get("stream", -1)),
    }
    print("streams:", streams)
    if sorted(streams.values()) != [1, 2, 3]:
        fail("expected disjoint streams 1,2,3", failures)

    print("\n== histogram (sampling) ==")
    bins = []
    with open(paths["histogram.csv"]) as f:
        for line in f:
            if line.startswith("#") or line.startswith("bin"):
                continue
            parts = line.strip().split(",")
            if len(parts) == 4:
                bins.append(int(parts[3]))
    declared = int(hist_meta["draws"])
    total = sum(bins)
    nbins = len(bins)
    expected = declared / nbins
    max_rel = max(abs(c - expected) / expected for c in bins)
    chi2 = sum((c - expected) ** 2 / expected for c in bins)
    print(f"bins={nbins} sum={total} declared_draws={declared} "
          f"chi2={chi2:.4f} max_rel={max_rel:.5f}")
    if total != declared:
        fail("histogram counts do not sum to declared draws", failures)
    if max_rel >= 0.05 or chi2 >= 40.0:
        fail("histogram fails uniformity acceptance", failures)

    print("\n== variance (material) ==")
    mean = var_meta["mean"]
    var = var_meta["variance"]
    analytic_var = var_meta["analytic_variance"]
    print(f"mean={mean:.9f} variance={var:.9f} "
          f"analytic_variance={analytic_var:.9f} "
          f"abs_err={abs(var - analytic_var):.2e}")
    if not var_meta.get("pass"):
        fail("variance artifact self-reports failure", failures)
    if abs(var - analytic_var) >= 0.02:
        fail("variance outside tolerance of analytic value", failures)

    print("\n== image hash (geometry) ==")
    actual = f"{fnv1a64(paths['image.ppm']):016x}"
    recorded = img_meta["fnv1a64"]
    print(f"recomputed={actual} recorded={recorded}")
    if actual != recorded:
        fail("image hash mismatch", failures)

    print("\n== summary ==")
    if failures:
        for m in failures:
            print(f"FAIL: {m}")
        print("RESULT: FAIL")
        return 1
    print("all four artifacts present, same scene, three streams, "
          "histogram/variance/hash checks pass")
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python3 -I verify_artifacts.py <artdir>", file=sys.stderr)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
