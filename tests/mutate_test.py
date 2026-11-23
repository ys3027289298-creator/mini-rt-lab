#!/usr/bin/env python3
"""
Mutation (fault-injection) harness. Proves the minimal test entries actually
catch defects. No third-party libraries (stdlib only).

For each mutation it copies the module source into an isolated sandbox,
applies one exact text replacement (a deliberately injected defect),
rebuilds the corresponding test entry against the mutated source, and runs
it. A mutation is KILLED when the test exits non-zero; it SURVIVES when the
test still passes (a coverage gap).

Usage: python3 -I mutate_test.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
CC = os.environ.get("CC", "gcc")
CFLAGS = ["-O2", "-std=c11", "-Wall", "-Wextra"]

# (id, source file, description, find, replace). Anchors are exact source
# substrings; each must occur exactly once so a bad match cannot silently
# hit the wrong line.
MUTATIONS = [
    # ---------------- geometry / intersection path (bvh.c) ----------------
    ("G1", "bvh/src/bvh.c", "aabb_hit: swap the slab t0/t1 update bounds",
     "            if (ta > t0) t0 = ta;\n            if (tb < t1) t1 = tb;",
     "            if (tb > t0) t0 = tb;\n            if (ta < t1) t1 = ta;"),
    ("G2", "bvh/src/bvh.c", "aabb_hit: drop the early-exit t0>t1 rejection",
     "            if (t0 > t1) return 0;",
     "            if (t0 > t1 && 0) return 0;"),
    ("G3", "bvh/src/bvh.c", "aabb_hit: parallel-ray outside test inverted",
     "            if (o[i] < b->min[i] || o[i] > b->max[i]) return 0;",
     "            if (o[i] >= b->min[i] && o[i] <= b->max[i]) return 0;"),
    ("G4", "bvh/src/bvh.c", "aabb_grow: min/max comparison direction reversed",
     "        if (o->min[i] < b->min[i]) b->min[i] = o->min[i];\n        if (o->max[i] > b->max[i]) b->max[i] = o->max[i];",
     "        if (o->min[i] > b->min[i]) b->min[i] = o->min[i];\n        if (o->max[i] < b->max[i]) b->max[i] = o->max[i];"),
    ("G5", "bvh/src/bvh.c", "aabb_area: surface-area factor 2 dropped",
     "    return 2.0f * (dx * dy + dy * dz + dz * dx);",
     "    return 1.0f * (dx * dy + dy * dz + dz * dx);"),
    ("G6", "bvh/src/bvh.c", "build: widest-axis choice picks wrong axis order",
     "    int axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2) : (w[1] > w[2] ? 1 : 2);",
     "    int axis = (w[0] < w[1]) ? (w[0] < w[2] ? 0 : 2) : (w[1] < w[2] ? 1 : 2);"),
    ("G7", "bvh/src/bvh.c", "build: centroid split condition uses <= instead of <",
     "        if (v < c[axis]) {",
     "        if (v <= c[axis]) {"),
    ("G8", "bvh/src/bvh.c", "build: degenerate midpoint fallback removed (forces one-sided split)",
     "    if (mid == first || mid == first + count) mid = first + count / 2;",
     "    if (0 && (mid == first || mid == first + count)) mid = first + count / 2;"),
    ("G9", "bvh/src/bvh.c", "aabb_hit: returned entry distance is t1 instead of t0",
     "    if (tmin) *tmin = t0;",
     "    if (tmin) *tmin = t1;"),
    ("G10", "bvh/src/bvh.c", "aabb_hit: negative-direction slab handling broken (no ta/tb swap)",
     "            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }",
     "            if (ta < tb) { float tmp = ta; ta = tb; tb = tmp; }"),
    ("G11", "bvh/src/bvh.c", "aabb_hit: parallel tolerance loosened 100x",
     "        if (fabsf(d[i]) < 1e-8f) {",
     "        if (fabsf(d[i]) < 1e-6f) {"),
    ("G12", "bvh/src/bvh.c", "aabb_area: negative-extent guard removed",
     "    if (dx < 0 || dy < 0 || dz < 0) return 0.0f;",
     "    if (dx < 0 && dy < 0 && dz < 0) return 0.0f;"),

    # ---------------- material evaluation path (material.c) ----------------
    ("M1", "material/src/material.c", "Lambert: divide by pi^2 instead of pi",
     "static float brdf_lambert(float albedo) { return albedo / 3.14159265f; }",
     "static float brdf_lambert(float albedo) { return albedo / (3.14159265f*3.14159265f); }"),
    ("M2", "material/src/material.c", "GGX: a2 computed as a instead of a*a",
     "    float a2 = a * a;",
     "    float a2 = a;"),
    ("M3", "material/src/material.c", "GGX: pi dropped from denominator",
     "    return a2 / (3.14159265f * d * d);",
     "    return a2 / (d * d);"),
    ("M4", "material/src/material.c", "GGX: denominator uses (a2+1) instead of (a2-1)",
     "    float d = nh * nh * (a2 - 1.0f) + 1.0f;",
     "    float d = nh * nh * (a2 + 1.0f) + 1.0f;"),
    ("M5", "material/src/material.c", "cosine_pdf: factor of pi removed",
     "static float cosine_pdf(float cos_theta) { return cos_theta / 3.14159265f; }",
     "static float cosine_pdf(float cos_theta) { return cos_theta; }"),
    ("M6", "material/src/material.c", "cosine_sample: radial r=u1 instead of sqrt(u1)",
     "    float r = sqrtf(u1);",
     "    float r = u1;"),
    ("M7", "material/src/material.c", "cosine_sample: phi uses u1 instead of u2",
     "    float phi = 2.0f * 3.14159265f * u2;",
     "    float phi = 2.0f * 3.14159265f * u1;"),
    ("M8", "material/src/material.c", "cosine_sample: z uses u1 instead of 1-u1",
     "    V3 v = { r * cosf(phi), r * sinf(phi), sqrtf(1.0f - u1 > 0 ? 1.0f - u1 : 0) };",
     "    V3 v = { r * cosf(phi), r * sinf(phi), sqrtf(u1 > 0 ? u1 : 0) };"),
    ("M9", "material/src/material.c", "cosine_sample: normalization removed",
     "    return norm3(v);",
     "    return v;"),
    ("M10", "material/src/material.c", "GGX: numerator a2 replaced by a (constant scaling)",
     "    return a2 / (3.14159265f * d * d);",
     "    return a / (3.14159265f * d * d);"),
    ("M11", "material/src/material.c", "Lambert: additive +0.1 bias on the BRDF",
     "static float brdf_lambert(float albedo) { return albedo / 3.14159265f; }",
     "static float brdf_lambert(float albedo) { return albedo / 3.14159265f + 0.1f; }"),
    ("M12", "material/src/material.c", "cosine_sample: r sqrt argument clamped wrong (u1+0.25)",
     "    float r = sqrtf(u1);",
     "    float r = sqrtf(u1 + 0.25f);"),
]

TEST_FOR = {
    "bvh/src/bvh.c": ("bvh_test.c", "BVH_SRC", "bvh_mut"),
    "material/src/material.c": ("material_test.c", "MATERIAL_SRC", "material_mut"),
}


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def main():
    sandbox = tempfile.mkdtemp(prefix="mini_rt_mutation_")
    results = []
    try:
        for mid, rel, desc, find, repl in MUTATIONS:
            src_path = os.path.join(REPO, rel)
            with open(src_path, "r") as f:
                original = f.read()
            if original.count(find) != 1:
                results.append((mid, rel, desc, "ANCHOR-ERROR", []))
                continue
            mutated = original.replace(find, repl)
            mut_file = os.path.join(sandbox, rel.replace("/", "_"))
            with open(mut_file, "w") as f:
                f.write(mutated)

            test_c, macro, bin_name = TEST_FOR[rel]
            bin_path = os.path.join(sandbox, bin_name)
            compile_cmd = [CC] + CFLAGS + [
                '-D%s="%s"' % (macro, mut_file),
                os.path.join(ROOT, test_c), "-o", bin_path, "-lm",
            ]
            cp = run(compile_cmd)
            if cp.returncode != 0:
                # compile failure counts as killed (test build rejects defect),
                # but record distinctly
                results.append((mid, rel, desc, "KILLED-COMPILE",
                               cp.stderr.strip().splitlines()[:3]))
                continue
            rp = run([bin_path])
            failed = [ln for ln in rp.stdout.splitlines()
                      if ln.endswith("=0")]
            state = "KILLED" if rp.returncode != 0 else "SURVIVED"
            results.append((mid, rel, desc, state, failed))
    finally:
        shutil.rmtree(sandbox, ignore_errors=True)

    killed = [r for r in results if r[3].startswith("KILLED")]
    survived = [r for r in results if r[3] == "SURVIVED"]
    anchor_err = [r for r in results if r[3] == "ANCHOR-ERROR"]
    total = len(results)

    print("=" * 78)
    print("MUTATION TEST RESULTS")
    print("=" * 78)
    path_group = {}
    for mid, rel, desc, state, failed in results:
        path_group.setdefault(rel, []).append((mid, state))
        mark = "KILLED " if state.startswith("KILLED") else "SURVIVE"
        print("[%s] %-4s %-24s %s" % (mark, mid, rel.split("/")[0], desc))
        if failed:
            for fn in failed[:6]:
                print("          failing check -> %s" % fn)
            if len(failed) > 6:
                print("          ... and %d more" % (len(failed) - 6))

    print("-" * 78)
    for rel, rows in path_group.items():
        k = sum(1 for _, s in rows if s.startswith("KILLED"))
        n = len(rows)
        print("%-24s killed %d/%d  capture %.1f%%"
              % (rel, k, n, 100.0 * k / n))
    rate = 100.0 * len(killed) / total if total else 0.0
    print("-" * 78)
    print("TOTAL %d mutations: %d killed, %d survived, %d anchor errors"
          % (total, len(killed), len(survived), len(anchor_err)))
    print("OVERALL CAPTURE RATE: %.1f%%" % rate)
    if survived:
        print("\nSURVIVING MUTATIONS (coverage gaps):")
        for mid, rel, desc, state, _ in survived:
            print("  %s (%s): %s" % (mid, rel, desc))
    if rate < 80.0:
        print("\nWARNING: capture rate below 80%%; missing cases listed above.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
