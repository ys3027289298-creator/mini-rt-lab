#!/usr/bin/env python3
"""Image-oracle mutation testing for bvh/src/bvh.c and material/src/material.c.

For every mutant in mutants.json:
  1. patch a sandbox copy of the module source,
  2. compile the module's own self-test and run it (unit-test baseline),
  3. compile the renderer (mut/render) against the patched source,
  4. render every scene with the same camera/RNG as the pristine reference,
  5. diff the 8-bit frames and classify the mutant as caught/survived.

"Caught" means the *image* differs visibly: at least MIN_PIX pixels differ by
at least VIS_DIFF in some RGB channel.  Bitwise-only differences are recorded
separately, because they can come from harmless float reassociation.

Stdlib only.  Usage:
  tools/run_mutations.py [--jobs N] [--width W] [--height H] [--spp S]
                         [--seed K] [--scenes grid,rough,occlude]
                         [--out results] [--only id,id]
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
RENDER = os.path.join(ROOT, "render")
CATALOGUE = os.path.join(ROOT, "mutants.json")

VIS_DIFF = 3      # 8-bit channel delta treated as visible
MIN_PIX = 8       # pixels that must exceed it


def load_catalogue():
    with open(CATALOGUE, "r", encoding="utf-8") as f:
        return json.load(f)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_frame(path):
    with open(path, "rb") as f:
        d = f.read()
    if d[:4] != b"RT1\x00":
        raise SystemExit("%s: not an RT1 frame" % path)
    w, h, spp, hv = struct.unpack("<IIII", d[4:20])
    px = d[20:]
    if len(px) != w * h * 3:
        raise SystemExit("%s: truncated frame" % path)
    return w, h, spp, hv, px


def diff_stats(a, b):
    """Compare two RGB8 byte strings."""
    if len(a) != len(b):
        return {"max": 255, "n_vis": 1 << 30, "rmse": 255.0, "bitwise": True}
    maxd = 0
    nvis = 0
    se = 0
    for i in range(0, len(a), 3):
        d0 = abs(a[i] - b[i])
        d1 = abs(a[i + 1] - b[i + 1])
        d2 = abs(a[i + 2] - b[i + 2])
        m = d0 if d0 >= d1 and d0 >= d2 else (d1 if d1 >= d2 else d2)
        if m > maxd:
            maxd = m
        if m >= VIS_DIFF:
            nvis += 1
        se += d0 * d0 + d1 * d1 + d2 * d2
    n = max(1, len(a) // 3)
    return {"max": maxd, "n_vis": nvis, "rmse": (se / (3.0 * n)) ** 0.5,
            "bitwise": a != b}


def write_diff_raw(path, a, b):
    """Canonical diff artifact: grayscale 4*max-channel-delta, clipped."""
    out = bytearray(len(a) // 3)
    for i in range(0, len(a), 3):
        m = max(abs(a[i] - b[i]), abs(a[i + 1] - b[i + 1]), abs(a[i + 2] - b[i + 2]))
        out[i // 3] = min(255, m * 16)
    with open(path, "wb") as f:
        f.write(bytes(out))
    return hashlib.sha256(bytes(out)).hexdigest()


def run(cmd, cwd=None, binary=False):
    """Run a command; stdout/stderr come back as str unless binary=True."""
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    if binary:
        return p.returncode, p.stdout, p.stderr.decode("utf-8", "replace")
    return (p.returncode, p.stdout.decode("utf-8", "replace"),
            p.stderr.decode("utf-8", "replace"))


def compile_renderer(build_dir, bvh_src, mat_src, tag):
    bvh_src = os.path.abspath(bvh_src)
    mat_src = os.path.abspath(mat_src)
    exe = os.path.join(build_dir, "rt_" + tag)
    cmd = [os.environ.get("CC", "gcc"), "-O2", "-std=c11", "-Wall", "-Wextra",
           "-I" + RENDER,
           '-DBVH_SRC="%s"' % bvh_src,
           '-DMAT_SRC="%s"' % mat_src,
           os.path.join(RENDER, "rt.c"),
           os.path.join(RENDER, "bvh_shim.c"),
           os.path.join(RENDER, "mat_shim.c"),
           "-o", exe, "-lm"]
    rc, out, err = run(cmd)
    return rc, err, exe


def render(exe, scene, w, h, spp, seed, out_path):
    # frames are binary; capture as bytes (a mutant may emit arbitrary output)
    rc, out, err = run([exe, scene, str(w), str(h), str(spp), str(seed)], binary=True)
    if rc != 0:
        return rc, err, None
    with open(out_path, "wb") as f:
        f.write(out)
    return 0, err, out_path


def apply_mutant(mut, src_text):
    find = mut["find"]
    if src_text.count(find) != 1:
        raise RuntimeError("find-string matched %d times" % src_text.count(find))
    return src_text.replace(find, mut["replace"], 1)


def draw_png(path, w, h, gray):
    import zlib
    raw = b"".join(b"\x00" + gray[y * w:(y + 1) * w] for y in range(h))

    def chunk(t, d):
        c = t + d
        return struct.pack(">I", len(d)) + c + struct.pack(">I", zlib.crc32(c) & 0xffffffff)

    png = b"\x89PNG\r\n\x1a\n"
    png += chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0))
    png += chunk(b"IDAT", zlib.compress(raw, 9))
    png += chunk(b"IEND", b"")
    with open(path, "wb") as f:
        f.write(png)


def process_mutant(mut, cat, cfg, refs):
    mid = mut["id"]
    sandbox = os.path.join(cfg["out"], "mutants", mid)
    os.makedirs(sandbox, exist_ok=True)
    root = cat["root"]
    bvh_src = os.path.join(root, "bvh", "src", "bvh.c")
    mat_src = os.path.join(root, "material", "src", "material.c")
    rec = {
        "id": mid, "module": mut["module"], "operator": mut["operator"],
        "kind": mut["kind"], "title": mut["title"],
        "equivalence_note": mut.get("equivalence_note"),
        "selftest": {"rc": None, "caught": None, "failing": []},
        "build": {"rc": None, "error": None},
        "scenes": {}, "caught": False, "caught_bitwise": False,
        "survived": None, "max_diff": 0, "max_diff_scene": None,
        "reason": None, "source_sha256": None, "sandbox": None,
        "expected": mut.get("expected", "caught"), "as_expected": None,
    }

    # --- unit-test baseline: compile the mutated module as its own program ---
    target = os.path.join(sandbox, os.path.basename(mut["module"] + ".c"))
    with open(bvh_src if mut["module"] == "bvh" else mat_src, "r", encoding="utf-8") as f:
        text = f.read()
    try:
        patched = apply_mutant(mut, text)
    except RuntimeError as e:
        rec["reason"] = "patch failed: %s" % e
        rec["as_expected"] = False
        return rec
    with open(target, "w", encoding="utf-8") as f:
        f.write(patched)
    rec["source_sha256"] = hashlib.sha256(patched.encode("utf-8")).hexdigest()
    rec["sandbox"] = os.path.relpath(sandbox, cfg["out"])
    st_exe = os.path.join(sandbox, "selftest")
    rc, _, err = run([os.environ.get("CC", "gcc"), "-O2", "-std=c11", "-Wall",
                      "-Wextra", target, "-o", st_exe, "-lm"])
    if rc != 0:
        rec["selftest"] = {"rc": rc, "caught": True, "failing": ["compile-error"]}
    else:
        rc2, out2, _ = run([st_exe])
        failing = [ln.split("=")[0] for ln in out2.splitlines()
                   if ln.endswith("=0")]
        rec["selftest"] = {"rc": rc2, "caught": rc2 != 0, "failing": failing}

    # --- renderer build ---
    bvh_use = target if mut["module"] == "bvh" else bvh_src
    mat_use = target if mut["module"] == "material" else mat_src
    rc, err, exe = compile_renderer(sandbox, bvh_use, mat_use, "mut")
    rec["build"] = {"rc": rc, "error": err.strip()[:2000] if rc else None}
    if rc != 0:
        rec["reason"] = "renderer build failed"
        rec["as_expected"] = False
        return rec

    # --- render + diff ---
    for scene in cfg["scenes"]:
        outp = os.path.join(sandbox, scene + ".raw")
        rc, err, _ = render(exe, scene, cfg["w"], cfg["h"], cfg["spp"],
                            cfg["seed"], outp)
        if rc != 0:
            rec["scenes"][scene] = {"status": "render_error", "detail": err.strip()[:500]}
            # a crash / poisoned frame on valid input is observable too
            rec["scenes"][scene]["diff"] = {"max": 255, "n_vis": 1 << 30,
                                            "rmse": 255.0, "bitwise": True}
            rec["caught"] = True
            continue
        w, h, spp, hv, px = read_frame(outp)
        ref = refs[scene]
        st = diff_stats(px, ref["px"])
        diff_raw = os.path.join(sandbox, "diff_%s.raw" % scene)
        dh = write_diff_raw(diff_raw, px, ref["px"])
        draw_png(os.path.join(sandbox, "diff_%s.png" % scene), w, h,
                 open(diff_raw, "rb").read())
        st.update({
            "status": "ok",
            "frame_sha256": hashlib.sha256(px).hexdigest(),
            "diff_sha256": dh,
            "diff_png": "mutants/%s/diff_%s.png" % (mid, scene),
            "ref_differs": st["bitwise"],
        })
        rec["scenes"][scene] = st
        if st["bitwise"]:
            rec["caught_bitwise"] = True
        if st["max"] >= VIS_DIFF and st["n_vis"] >= MIN_PIX:
            rec["caught"] = True
        if st["max"] > rec["max_diff"]:
            rec["max_diff"] = st["max"]
            rec["max_diff_scene"] = scene

    expect = mut.get("expected", "caught")
    rec["expected"] = expect
    rec["survived"] = not rec["caught"]
    rec["as_expected"] = (rec["caught"] and expect == "caught") or \
                         (not rec["caught"] and expect == "survive")
    if rec["kind"] == "equivalent":
        rec["reason"] = ("equivalent mutant: visually identical, as expected"
                         if not rec["caught"] else
                         "WRONG: equivalent mutant produced a visible diff")
    elif rec["kind"] == "defect-subthreshold":
        rec["reason"] = ("below the measured sensitivity floor, as expected"
                         if not rec["caught"] else
                         "caught despite being predicted sub-threshold")
    elif rec["caught"]:
        rec["reason"] = "caught by image oracle"
    else:
        rec["reason"] = "UNEXPECTED SURVIVOR: defect invisible to all scenes"
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=240)
    ap.add_argument("--spp", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scenes", default="grid,rough,occlude,grazing,warm")
    ap.add_argument("--out", default=os.path.join(ROOT, "results"))
    ap.add_argument("--only", default="")
    ap.add_argument("--build", default=os.path.join(ROOT, "build"))
    a = ap.parse_args()

    cat = load_catalogue()
    root = cat["root"]
    cfg = {"w": a.width, "h": a.height, "spp": a.spp, "seed": a.seed,
           "scenes": [s for s in a.scenes.split(",") if s], "build": a.build,
           "out": a.out}
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(os.path.join(a.out, "frames"), exist_ok=True)
    os.makedirs(os.path.join(a.build, "ref"), exist_ok=True)

    bvh_src = os.path.join(root, "bvh", "src", "bvh.c")
    mat_src = os.path.join(root, "material", "src", "material.c")

    # reference renderer + frames
    rc, err, ref_exe = compile_renderer(os.path.join(a.build, "ref"), bvh_src, mat_src, "ref")
    if rc != 0:
        sys.stderr.write("reference build failed:\n" + err)
        return 1
    refs = {}
    ref_meta = {}
    for scene in cfg["scenes"]:
        outp = os.path.join(a.build, "ref", scene + ".raw")
        rc, err, _ = render(ref_exe, scene, a.width, a.height, a.spp, a.seed, outp)
        if rc != 0:
            sys.stderr.write("reference render failed for %s:\n%s\n" % (scene, err))
            return 1
        w, h, spp, hv, px = read_frame(outp)
        refs[scene] = {"px": px}
        ref_meta[scene] = {
            "sha256": hashlib.sha256(px).hexdigest(),
            "fnv1a_header": "%08x" % hv,
            "raw": "frames/%s.raw" % scene,
            "bytes": len(px),
        }
        shutil.copyfile(outp, os.path.join(a.out, "frames", scene + ".raw"))

    mutants = cat["mutants"]
    if a.only:
        want = set(a.only.split(","))
        mutants = [m for m in mutants if m["id"] in want]

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(process_mutant, m, cat, cfg, refs): m for m in mutants}
        for i, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            rec = fut.result()
            results.append(rec)
            mark = "CAUGHT" if rec["caught"] else ("SURVIVED" if rec["kind"] == "defect" else "survived(equiv)")
            print("[%d/%d] %-28s %-9s %s" % (i, len(mutants), rec["id"], rec["kind"], mark),
                  flush=True)

    order = {m["id"]: i for i, m in enumerate(mutants)}
    results.sort(key=lambda r: order.get(r["id"], 1 << 30))

    def rate(kind):
        rows = [r for r in results if r["kind"] == kind]
        caught = sum(1 for r in rows if r["caught"])
        return {"n": len(rows), "caught": caught,
                "rate": (caught / len(rows)) if rows else None}

    defects = [r for r in results if r["kind"] in ("defect", "defect-subthreshold")]
    summary = {
        "defects": rate("defect"),
        "defects_incl_subthreshold": {
            "n": len(defects),
            "caught": sum(1 for r in defects if r["caught"]),
            "rate": (sum(1 for r in defects if r["caught"]) / len(defects)) if defects else None,
        },
        "defect_subthreshold": rate("defect-subthreshold"),
        "equivalents": rate("equivalent"),
        "as_expected": sum(1 for r in results if r["as_expected"]),
        "wrong": [r["id"] for r in results if not r["as_expected"]],
        "by_operator": {},
        "by_module": {},
    }
    for key, field in (("by_operator", "operator"), ("by_module", "module")):
        for r in results:
            b = summary[key].setdefault(r[field], {"n": 0, "caught": 0})
            if r["kind"] in ("defect", "defect-subthreshold"):
                b["n"] += 1
                b["caught"] += 1 if r["caught"] else 0
    out = {
        "config": cfg,
        "thresholds": {"vis_diff": VIS_DIFF, "min_pixels": MIN_PIX},
        "reference": ref_meta,
        "summary": summary,
        "results": results,
        "compiler": run([os.environ.get("CC", "gcc"), "--version"])[1].splitlines()[0],
    }
    with open(os.path.join(a.out, "results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print("wrote", os.path.join(a.out, "results.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
