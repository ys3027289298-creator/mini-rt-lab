#!/usr/bin/env python3
"""Diff rendered frames between a git revision and the current working tree.

Answers "did my edits to bvh/material change the picture?" without guessing:
both revisions are compiled into the same renderer and compared on identical
cameras, resolutions, sample counts and seeds.  Because rendering is paired and
deterministic (see tools/validate.py), any nonzero visible difference is real.

Usage:
  tools/compare_rev.py [--rev HEAD] [--width 320] [--height 240] [--spp 16]
                       [--seed 1] [--scenes grid,rough,occlude,grazing,warm]

Exit 0 if all scenes are identical, 1 if any differ, 2 on error.
Stdlib only.
"""
import argparse
import hashlib
import json
import os
import struct
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
REPO = os.path.abspath(os.path.join(ROOT, os.pardir))
RENDER = os.path.join(ROOT, "render")
VIS_DIFF = 3
MIN_PIX = 8


def run(cmd, env=None, binary=False, cwd=None):
    e = dict(os.environ)
    if env:
        e.update(env)
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, env=e)
    out = p.stdout if binary else p.stdout.decode("utf-8", "replace")
    return p.returncode, out, p.stderr.decode("utf-8", "replace")


def git_show(rev, path, dest):
    rc, data, err = run(["git", "show", "%s:%s" % (rev, path)], cwd=REPO, binary=True)
    if rc != 0:
        raise SystemExit("git show %s:%s failed:\n%s" % (rev, path, err))
    with open(dest, "wb") as f:
        f.write(data)
    return dest


def build(tag, bvh_src, mat_src):
    d = os.path.join(ROOT, "build", "compare", tag)
    os.makedirs(d, exist_ok=True)
    exe = os.path.join(d, "rt")
    rc, _, err = run([os.environ.get("CC", "gcc"), "-O2", "-std=c11", "-Wall",
                      "-Wextra", "-I" + RENDER,
                      '-DBVH_SRC="%s"' % os.path.abspath(bvh_src),
                      '-DMAT_SRC="%s"' % os.path.abspath(mat_src),
                      os.path.join(RENDER, "rt.c"),
                      os.path.join(RENDER, "bvh_shim.c"),
                      os.path.join(RENDER, "mat_shim.c"),
                      "-o", exe, "-lm"])
    if rc != 0:
        sys.stderr.write(err)
        raise SystemExit("build failed for " + tag)
    return exe


def frame(p):
    d = open(p, "rb").read()
    w, h, spp, hv = struct.unpack("<IIII", d[4:20])
    return d[20:]


def cmp(a, b):
    md = nv = 0
    for i in range(0, len(a), 3):
        m = max(abs(a[i] - b[i]), abs(a[i + 1] - b[i + 1]), abs(a[i + 2] - b[i + 2]))
        md = max(md, m)
        if m >= VIS_DIFF:
            nv += 1
    return md, nv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rev", default="HEAD")
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=240)
    ap.add_argument("--spp", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scenes", default="grid,rough,occlude,grazing,warm")
    a = ap.parse_args()
    scenes = [s for s in a.scenes.split(",") if s]

    rc, head, _ = run(["git", "rev-parse", "--short", a.rev], cwd=REPO)
    if rc != 0:
        raise SystemExit("bad revision: " + a.rev)
    head = head.strip()

    with tempfile.TemporaryDirectory() as td:
        bvh_old = git_show(a.rev, "bvh/src/bvh.c", os.path.join(td, "bvh.c"))
        mat_old = git_show(a.rev, "material/src/material.c", os.path.join(td, "material.c"))
        old = build("rev_" + head, bvh_old, mat_old)
    new = build("worktree",
                os.path.join(REPO, "bvh", "src", "bvh.c"),
                os.path.join(REPO, "material", "src", "material.c"))

    outdir = os.path.join(ROOT, "results", "compare")
    os.makedirs(outdir, exist_ok=True)
    report = {"rev": head, "config": vars(a), "scenes": {}}
    differs = False
    for s in scenes:
        args = [s, str(a.width), str(a.height), str(a.spp), str(a.seed)]
        _, po, _ = run([old, *args], binary=True)
        _, pn, _ = run([new, *args], binary=True)
        fo = os.path.join(td, "old.raw")
        fn = os.path.join(outdir, s + ".raw")
        open(os.path.join(outdir, "rev_" + s + ".raw"), "wb").write(po)
        open(fn, "wb").write(pn)
        A, B = frame(os.path.join(outdir, "rev_" + s + ".raw")), frame(fn)
        md, nv = cmp(A, B)
        rec = {
            "visible_pixels": nv,
            "max_channel_delta": md,
            "identical_bytes": A == B,
            "rev_frame_sha256": hashlib.sha256(A).hexdigest(),
            "worktree_frame_sha256": hashlib.sha256(B).hexdigest(),
            "visibly_different": nv >= MIN_PIX and md >= VIS_DIFF,
        }
        diff = bytearray(len(A) // 3)
        for i in range(0, len(A), 3):
            m = max(abs(A[i] - B[i]), abs(A[i + 1] - B[i + 1]), abs(A[i + 2] - B[i + 2]))
            diff[i // 3] = min(255, m * 16)
        dp = os.path.join(outdir, "diff_" + s + ".raw")
        open(dp, "wb").write(bytes(diff))
        rec["diff_sha256"] = hashlib.sha256(bytes(diff)).hexdigest()
        report["scenes"][s] = rec
        if rec["visibly_different"]:
            differs = True
        print("%-8s %s: visible_px=%d max_delta=%d%s" %
              (s, "CHANGED" if rec["visibly_different"] else "same", nv, md,
               "" if A != B else " (byte-identical)"),
              flush=True)

    json.dump(report, open(os.path.join(outdir, "compare.json"), "w",
                           encoding="utf-8"), indent=2)
    print("wrote", os.path.join(outdir, "compare.json"))
    return 1 if differs else 0


if __name__ == "__main__":
    sys.exit(main())
