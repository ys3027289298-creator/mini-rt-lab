#!/usr/bin/env python3
"""Measure the oracle's sensitivity floor.

Takes the Lambert denominator (the cleanest scalar in the material module) and
perturbs pi by a geometric series of relative errors, rendering all scenes each
time.  Reports the smallest relative error that still produces a visible frame
difference, and the error at which detection becomes robust across every scene.

This bounds the claim "the image oracle catches material defects": it can see
scaled errors above the floor, and is blind below it.  The floor is a property
of 8-bit quantization plus scene exposure, not of the mutation harness.

Stdlib only.  Writes sensitivity.json.
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
VIS_DIFF = 3
MIN_PIX = 8
PI = "3.14159265f"


def run(cmd, env=None, binary=False):
    e = dict(os.environ)
    if env:
        e.update(env)
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=e)
    return p.returncode, (p.stdout if binary else p.stdout.decode("utf-8", "replace")), \
        p.stderr.decode("utf-8", "replace")


def frame_px(blob):
    return blob[20:]


def vis(a, b):
    if len(a) != len(b):
        return 1 << 30, 255
    nv = 0
    md = 0
    for i in range(0, len(a), 3):
        m = max(abs(a[i] - b[i]), abs(a[i + 1] - b[i + 1]), abs(a[i + 2] - b[i + 2]))
        md = max(md, m)
        if m >= VIS_DIFF:
            nv += 1
    return nv, md


def build(tag, bvh_src, mat_src):
    d = os.path.join(ROOT, "build", "sensitivity", tag)
    os.makedirs(d, exist_ok=True)
    exe = os.path.join(d, "rt")
    rc, _, err = run([os.environ.get("CC", "gcc"), "-O2", "-std=c11", "-Wall",
                      "-Wextra", "-I" + RENDER,
                      '-DBVH_SRC="%s"' % bvh_src,
                      '-DMAT_SRC="%s"' % mat_src,
                      os.path.join(RENDER, "rt.c"),
                      os.path.join(RENDER, "bvh_shim.c"),
                      os.path.join(RENDER, "mat_shim.c"),
                      "-o", exe, "-lm"])
    if rc != 0:
        sys.stderr.write(err)
        raise SystemExit("build failed: " + tag)
    return exe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=240)
    ap.add_argument("--spp", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scenes", default="grid,rough,occlude,grazing,warm")
    ap.add_argument("--epsilons",
                    default="0.0001,0.0003,0.001,0.002,0.003,0.005,0.01,0.02,0.05,0.1")
    a = ap.parse_args()
    cat = json.load(open(os.path.join(ROOT, "mutants.json"), encoding="utf-8"))
    root = cat["root"]
    bvh = os.path.join(root, "bvh", "src", "bvh.c")
    mat = os.path.join(root, "material", "src", "material.c")
    scenes = [s for s in a.scenes.split(",") if s]
    eps = [float(x) for x in a.epsilons.split(",")]

    base = open(mat, encoding="utf-8").read()
    anchor = "static float brdf_lambert(float albedo) { return albedo / %s; }" % PI
    if base.count(anchor) != 1:
        raise SystemExit("lambert anchor not found exactly once")

    ref_exe = build("ref", bvh, mat)
    refs = {}
    for s in scenes:
        _, out, _ = run([ref_exe, s, str(a.width), str(a.height), str(a.spp),
                         str(a.seed)], binary=True)
        refs[s] = frame_px(out)

    rows = []
    for e in eps:
        # albedo / pi  ->  albedo / (pi*(1+e))  ==  multiply the BRDF by 1/(1+e).
        # Keep the code shape identical to the pristine line otherwise.
        patched = base.replace(anchor, "static float brdf_lambert(float albedo) "
                                       "{ return albedo / (%s * %.9ff); }" % (PI, 1.0 + e), 1)
        p = os.path.join(ROOT, "build", "sensitivity", "src_%.6f.c" % e)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, "w", encoding="utf-8").write(patched)
        exe = build("e%.6f" % e, bvh, p)
        row = {"rel_error": e, "scenes": {}, "min_vis": None, "all_scenes_visible": True}
        for s in scenes:
            rc, out, err = run([exe, s, str(a.width), str(a.height), str(a.spp),
                                str(a.seed)], binary=True)
            if rc != 0:
                row["scenes"][s] = {"status": "error", "detail": err[:300]}
                continue
            nv, md = vis(frame_px(out), refs[s])
            row["scenes"][s] = {"n_vis": nv, "max": md,
                                "visible": nv >= MIN_PIX and md >= VIS_DIFF}
            if not row["scenes"][s]["visible"]:
                row["all_scenes_visible"] = False
        nvs = [v["n_vis"] for v in row["scenes"].values() if "n_vis" in v]
        row["min_vis"] = min(nvs) if nvs else None
        row["max_vis"] = max(nvs) if nvs else None
        rows.append(row)
        print("eps=%-8g min_vis=%-7s all_scenes_visible=%s" %
              (e, row["min_vis"], row["all_scenes_visible"]), flush=True)

    first_any = next((r for r in rows if any(v.get("visible") for v in r["scenes"].values())), None)
    first_all = next((r for r in rows if r["all_scenes_visible"]), None)
    out = {
        "probe": "brdf_lambert pi scaled by (1+eps); relative BRDF error = 1-1/(1+eps)",
        "config": vars(a),
        "thresholds": {"vis_diff": VIS_DIFF, "min_pixels": MIN_PIX},
        "rows": rows,
        "smallest_visible_any_scene": first_any["rel_error"] if first_any else None,
        "smallest_visible_in_every_scene": first_all["rel_error"] if first_all else None,
    }
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, "results", "sensitivity.json"), "w",
                        encoding="utf-8"), indent=2)
    print("smallest visible in any scene:", out["smallest_visible_any_scene"])
    print("smallest visible in every scene:", out["smallest_visible_in_every_scene"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
