#!/usr/bin/env python3
"""End-to-end consistency check.

For each scene, render with both implementations at several scales,
classify every diff pixel, check hand-labeled recall (zero misses), and
assert that conclusions are stable when the scale changes.

Run from consistency/ (or via `make test`):
    PYTHONPATH="$(python3 -c 'import site; print(site.getusersitepackages())')" \
        python3 -E tools/run_check.py
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

from consist_check import (  # noqa: E402
    BOUNDARY, CLASSES, COLOR_EPS, DEPTH_REL_EPS, OCCLUSION, PRECISION,
    UNCLASSIFIED, classify, read_fbuf,
)
from consist_check.labels import recall_against_labels  # noqa: E402

# (scene, label file or None, resolutions). Scales share the same field of
# view, so geometry occupies the same normalized region at every scale.
JOBS = [
    ("clean", None, (160, 320, 640)),
    ("lab", "lab_labels.json", (160, 320, 640)),
]
BASE_SCALE = 320
# Count growth is geometric: filled regions ~ O(s^2), edge/band effects
# ~ O(s). A class present in the scene must keep a scaling exponent in this
# physically bounded band at every scale pair, or the rule is flaky.
EXPONENT_MIN, EXPONENT_MAX = 0.5, 2.5
# The largest class (by count) must be the same at every scale.
EXAMPLE_PER_CLASS = 5


def render(impl, scene, w, h, out_path):
    exe = ROOT / impl
    scene_path = ROOT / "scenes" / f"{scene}.scene"
    subprocess.run(
        [str(exe), str(scene_path), str(w), str(h), str(out_path)],
        check=True, stdout=subprocess.DEVNULL,
    )


def examples_for(labels, metrics, rast, trace, cls, limit):
    ys, xs = np.nonzero(labels == cls)
    if len(ys) == 0:
        return []
    # rank by how strongly the pixel exhibits the class signature, so the
    # examples stay representative rather than landing on borderline cases
    if cls == PRECISION:
        score = np.maximum(metrics["color_delta"], metrics["depth_delta"])
    else:
        score = metrics["color_delta"]
    sel = np.argsort(score[ys, xs])[::-1][:limit]
    out = []
    for i in sel:
        y, x = int(ys[i]), int(xs[i])
        out.append({
            "x": x, "y": y,
            "id_rast": int(rast[y, x, 4]),
            "id_trace": int(trace[y, x, 4]),
            "color_delta": round(float(metrics["color_delta"][y, x]), 6),
            "depth_rel_delta": round(float(metrics["depth_delta"][y, x]), 6),
        })
    return out


def evaluate_scene(scene, label_name, sizes, out_dir):
    per_scale = {}
    for w in sizes:
        h = int(round(w * 0.75))
        rast_path = out_dir / f"{scene}_{w}x{h}_rast.fbuf"
        trace_path = out_dir / f"{scene}_{w}x{h}_trace.fbuf"
        render("raster", scene, w, h, rast_path)
        render("tracer", scene, w, h, trace_path)

        rast = read_fbuf(rast_path)
        trace = read_fbuf(trace_path)
        labels, metrics = classify(rast, trace)
        total_px = w * h
        counts = {cls: int(np.count_nonzero(labels == cls)) for cls in CLASSES}
        diff_total = sum(counts.values())
        per_scale[w] = {
            "height": h,
            "total_pixels": total_px,
            "diff_pixels": diff_total,
            "diff_fraction": diff_total / total_px,
            "counts": counts,
            "shares": {
                cls: (counts[cls] / diff_total if diff_total else 0.0)
                for cls in CLASSES
            },
            "labels": labels,
            "metrics": metrics,
            "rast": rast,
            "trace": trace,
        }

        if counts[UNCLASSIFIED] != 0:
            raise AssertionError(
                f"{scene} {w}x{h}: {counts[UNCLASSIFIED]} unclassified diff pixels; "
                "the taxonomy does not cover this disagreement"
            )

        save_diff_overlay(scene, w, h, labels, rast, out_dir)

    # hand-labeled recall gate at every scale
    recall = None
    if label_name is not None:
        recall = {}
        label_path = ROOT / "labels" / label_name
        for w in sizes:
            h = per_scale[w]["height"]
            seen, misses = recall_against_labels(
                out_dir / f"{scene}_{w}x{h}_rast.fbuf",
                out_dir / f"{scene}_{w}x{h}_trace.fbuf",
                label_path,
            )
            recall[w] = {"per_class": seen, "misses": misses,
                         "miss_count": len(misses)}
            if misses:
                raise AssertionError(
                    f"{scene} {w}x{h}: {len(misses)} hand-labeled pixels "
                    f"misclassified, first: {misses[0]}"
                )

    return per_scale, recall


def save_diff_overlay(scene, w, h, labels, rast, out_dir):
    color = {
        OCCLUSION: (231, 76, 60),    # red
        PRECISION: (46, 204, 113),   # green
        BOUNDARY: (241, 196, 15),    # yellow
        UNCLASSIFIED: (155, 89, 182),  # purple
    }
    img = (rast[:, :, 0:3] * 255.0).clip(0, 255).astype(np.uint8)
    img = np.repeat(img, 1, axis=0)
    for cls, rgb in color.items():
        mask = labels == cls
        img[mask] = rgb
    Image.fromarray(np.ascontiguousarray(img)).save(
        out_dir / f"{scene}_{w}x{h}_diff.png"
    )


def stability_assertions(scene, per_scale):
    """The conclusion must be invariant to resolution.

    Three gates, all derived from how pixel counts geometrically scale:
      1. the set of observed classes is identical at every scale;
      2. for each class, the count grows ~ s**k with k in the physically
         bounded band [EXPONENT_MIN, EXPONENT_MAX] (edges ~ O(s), areas
         ~ O(s^2)); a class that appears/disappears or explodes fails;
      3. the dominant class by count is the same at every scale.
    Per-class shares are reported but not gated, because an O(s) class
    naturally loses share to O(s^2) classes as resolution grows.
    """
    sizes = sorted(per_scale)
    present = {cls for cls in CLASSES if cls != UNCLASSIFIED
               and per_scale[BASE_SCALE]["counts"][cls] > 0}

    for w in sizes:
        cur_present = {cls for cls in CLASSES if cls != UNCLASSIFIED
                       and per_scale[w]["counts"][cls] > 0}
        if cur_present != present:
            raise AssertionError(
                f"{scene}: class set changed between scales "
                f"{sorted(present)} -> {sorted(cur_present)} at {w}px"
            )

    exponents = {}
    sw = np.array(sizes, dtype=float)
    for cls in sorted(present):
        counts = np.array([per_scale[w]["counts"][cls] for w in sizes], dtype=float)
        k = float(np.polyfit(np.log(sw), np.log(counts), 1)[0])
        exponents[cls] = round(k, 3)
        if not (EXPONENT_MIN <= k <= EXPONENT_MAX):
            raise AssertionError(
                f"{scene}: scaling exponent {k:.2f} for {cls} outside "
                f"[{EXPONENT_MIN}, {EXPONENT_MAX}] — classification is "
                "not stable across resolutions"
            )

    def dominant(w):
        return max(present, key=lambda c: per_scale[w]["counts"][c])
    leaders = {dominant(w) for w in sizes}
    if len(leaders) != 1:
        raise AssertionError(
            f"{scene}: dominant class changed across scales: {leaders}"
        )
    return exponents


def main():
    out_dir = ROOT / "out"
    out_dir.mkdir(exist_ok=True)
    report = {"scenes": {}, "thresholds": {
        "color_eps": COLOR_EPS, "depth_rel_eps": DEPTH_REL_EPS,
    }}
    failures = []

    for scene, label_name, sizes in JOBS:
        try:
            per_scale, recall = evaluate_scene(scene, label_name, sizes, out_dir)
            exponents = stability_assertions(scene, per_scale)
        except (AssertionError, subprocess.CalledProcessError) as exc:
            failures.append(f"{scene}: {exc}")
            continue

        base = per_scale[BASE_SCALE]
        entry = {
            "scales": {},
            "scaling_exponents": exponents,
            "examples_at_base_scale": {
                cls: examples_for(
                    base["labels"], base["metrics"],
                    base["rast"], base["trace"], cls, EXAMPLE_PER_CLASS,
                )
                for cls in CLASSES if base["counts"][cls] > 0
            },
        }
        for w in sizes:
            cur = per_scale[w]
            entry["scales"][f"{w}x{cur['height']}"] = {
                "diff_pixels": cur["diff_pixels"],
                "total_pixels": cur["total_pixels"],
                "diff_fraction": round(cur["diff_fraction"], 8),
                "counts": cur["counts"],
                "shares": {k: round(v, 6) for k, v in cur["shares"].items()},
            }
        if recall is not None:
            entry["label_recall"] = {
                str(w): {
                    "miss_count": recall[w]["miss_count"],
                    "per_class": recall[w]["per_class"],
                } for w in sizes
            }
        report["scenes"][scene] = entry

    report_path = out_dir / "report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    print_human_report(report, failures)
    if failures:
        print("\nRESULT: FAIL")
        for msg in failures:
            print(f"  - {msg}")
        return 1
    print("\nRESULT: PASS (zero unclassified, zero missed hand labels, "
          "conclusions stable across scales)")
    return 0


def print_human_report(report, failures):
    for name, entry in report["scenes"].items():
        print(f"\n=== scene: {name} ===")
        for scale, stats in entry["scales"].items():
            print(f"  {scale}: diff={stats['diff_pixels']} "
                  f"({stats['diff_fraction']:.6f} of pixels)")
            parts = [f"{cls}={stats['counts'][cls]}"
                     for cls in CLASSES if stats["counts"][cls] > 0]
            print(f"    by class: {', '.join(parts)}")
        print(f"  scaling exponents: {entry['scaling_exponents']}")
        if "label_recall" in entry:
            for scale, rec in entry["label_recall"].items():
                gate = "OK" if rec["miss_count"] == 0 else "MISS"
                detail = ", ".join(
                    f"{cls} {v['hit']}/{v['labeled']}"
                    for cls, v in rec["per_class"].items() if v["labeled"]
                )
                print(f"  label recall {scale}: {gate} ({detail})")
        for cls, exs in entry["examples_at_base_scale"].items():
            if not exs:
                continue
            coords = ", ".join(f"({e['x']},{e['y']})" for e in exs)
            print(f"  examples {cls}: {coords}")
    if failures:
        print("\nFailures:")
        for msg in failures:
            print(f"  - {msg}")


if __name__ == "__main__":
    sys.exit(main())
