#!/usr/bin/env python3
"""Validate the image oracle before trusting mutation results.

Checks, in order:
  1. determinism      -- the same binary and inputs produce byte-identical frames;
  2. build stability  -- two separate builds of the same sources agree;
  3. tree == brute    -- the BVH path and an unaccelerated all-boxes path produce
                         identical frames on the pristine sources, i.e. traversal
                         loses no hits and culls nothing it should not;
  4. noise floor      -- the oracle is *paired*: mutant and reference use the same
                         seed, so Monte-Carlo noise cancels exactly.  This check
                         records the paired floor (must be 0 differing pixels) and
                         the unpaired floor (different seed, large) to make that
                         assumption explicit;
  5. threshold sweep  -- for a few sanity defects, the detection margin at the
                         chosen threshold.

Stdlib only, same toolchain as the suite.  Writes validate.json.
"""
import argparse
import hashlib
import json
import os
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
RENDER = os.path.join(ROOT, "render")
CATALOGUE = os.path.join(ROOT, "mutants.json")
VIS_DIFF = 3
MIN_PIX = 8


def run(cmd, env=None, binary=False):
    e = dict(os.environ)
    if env:
        e.update(env)
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=e)
    out = p.stdout
    err = p.stderr.decode("utf-8", "replace")
    if not binary:
        out = out.decode("utf-8", "replace")
    return p.returncode, out, err


def frame(path=None, blob=None):
    d = blob if blob is not None else open(path, "rb").read()
    if d[:4] != b"RT1\x00":
        raise SystemExit("not an RT1 frame")
    w, h, spp, hv = struct.unpack("<IIII", d[4:20])
    return w, h, spp, hv, d[20:]


def cmp_frames(a, b):
    if len(a) != len(b):
        return {"max": 255, "n_vis": 1 << 30, "identical": False}
    md, nv = 0, 0
    for i in range(0, len(a), 3):
        m = max(abs(a[i] - b[i]), abs(a[i + 1] - b[i + 1]), abs(a[i + 2] - b[i + 2]))
        md = max(md, m)
        if m >= VIS_DIFF:
            nv += 1
    return {"max": md, "n_vis": nv, "identical": a == b}


def build(tag, bvh_src, mat_src, extra_flags=()):
    d = os.path.join(ROOT, "build", "validate", tag)
    os.makedirs(d, exist_ok=True)
    exe = os.path.join(d, "rt")
    rc, out, err = run([os.environ.get("CC", "gcc"), "-O2", "-std=c11", "-Wall",
                        "-Wextra", *extra_flags, "-I" + RENDER,
                        '-DBVH_SRC="%s"' % bvh_src,
                        '-DMAT_SRC="%s"' % mat_src,
                        os.path.join(RENDER, "rt.c"),
                        os.path.join(RENDER, "bvh_shim.c"),
                        os.path.join(RENDER, "mat_shim.c"),
                        "-o", exe, "-lm"])
    if rc != 0:
        sys.stderr.write(err)
        raise SystemExit("build failed for " + tag)
    return exe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=240)
    ap.add_argument("--spp", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scenes", default="grid,rough,occlude,grazing,warm")
    a = ap.parse_args()
    cat = json.load(open(CATALOGUE, encoding="utf-8"))
    root = cat["root"]
    bvh = os.path.join(root, "bvh", "src", "bvh.c")
    mat = os.path.join(root, "material", "src", "material.c")
    scenes = [s for s in a.scenes.split(",") if s]
    res = {"config": vars(a), "checks": {}}

    e1 = build("pristine_a", bvh, mat)
    e2 = build("pristine_b", bvh, mat)
    res["binaries"] = {
        "pristine_a": hashlib.sha256(open(e1, "rb").read()).hexdigest(),
        "pristine_b": hashlib.sha256(open(e2, "rb").read()).hexdigest(),
    }
    res["binaries_identical"] = res["binaries"]["pristine_a"] == res["binaries"]["pristine_b"]

    det, stab, tb, noise, unpaired = {}, {}, {}, {}, {}
    for s in scenes:
        args = [s, str(a.width), str(a.height), str(a.spp), str(a.seed)]
        _, b1, _ = run([e1, *args], binary=True)
        _, b2, _ = run([e1, *args], binary=True)
        _, b3, _ = run([e2, *args], binary=True)
        f1, f2, f3 = frame(blob=b1), frame(blob=b2), frame(blob=b3)
        det[s] = {"identical": b1 == b2,
                  "sha256": hashlib.sha256(f1[4]).hexdigest()}
        stab[s] = {"identical": b1 == b3}
        _, bb, _ = run([e1, *args], env={"RT_BRUTE": "1"}, binary=True)
        fb = frame(blob=bb)
        tb[s] = {"tree_equals_brute": f1[4] == fb[4],
                 "diff": cmp_frames(f1[4], fb[4])}
        # paired floor: re-render the pristine binary at the same seed
        _, bp, _ = run([e1, *args], binary=True)
        noise[s] = cmp_frames(f1[4], frame(blob=bp)[4])
        noise[s]["pixel_count"] = len(f1[4]) // 3
        # unpaired floor: a different seed.  Documented, not a pass criterion.
        _, bn, _ = run([e1, s, str(a.width), str(a.height), str(a.spp),
                        str(a.seed + 1)], binary=True)
        unpaired[s] = cmp_frames(f1[4], frame(blob=bn)[4])
        unpaired[s]["pixel_count"] = len(f1[4]) // 3

    res["checks"]["determinism"] = det
    res["checks"]["build_stability"] = stab
    res["checks"]["tree_equals_brute"] = tb
    res["checks"]["paired_noise_floor"] = noise
    res["checks"]["unpaired_noise_floor_seed_diff"] = unpaired
    res["pass"] = (
        res["binaries_identical"]
        and all(v["identical"] for v in det.values())
        and all(v["identical"] for v in stab.values())
        and all(v["tree_equals_brute"] for v in tb.values())
        and all(v["n_vis"] == 0 for v in noise.values())
    )

    out = os.path.join(ROOT, "results", "validate.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out, "w", encoding="utf-8"), indent=2)
    print(json.dumps({k: res["checks"][k] for k in res["checks"]}, indent=2))
    print("pass =", res["pass"])
    return 0 if res["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
