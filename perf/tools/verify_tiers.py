#!/usr/bin/env python3
"""
Two-tier performance-budget acceptance for the three rendering paths.

Reads the artifacts written by `perf run <tier>` and proves, for BOTH tiers:

  * completeness   -- metrics.json + all three PPMs exist (a missing tier or a
                      missing dataset is a failure, not a warning)
  * work budget    -- scheduled == width*height*spp, measured == scheduled,
                      and both are <= the tier's work bound
  * time budget    -- every per-path median AND the tier total are <= bounds
  * memory budget  -- peak RSS <= the tier's bound
  * image hash     -- FNV-1a64 recomputed over the exact PPM bytes equals the
                      hash recorded in metrics.json
  * no-check mode  -- optional: PPMs produced with --no-budget-check are
                      byte-identical to the checked run's PPMs

On success it also writes REPORT.md with the per-tier timing table, peak
memory, image hashes, bounds and the exact reproduce commands.

Runs under `python3 -I`; stdlib only.
Exit: 0 all tiers complete and within budget, 1 any failure, 2 usage/IO.
"""
import json
import os
import sys

TIERS = ("low", "high")
PATHS = ("bvh", "sample", "material")


def fnv1a64(path):
    h = 14695981039346656037
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            for byte in chunk:
                h ^= byte
                h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h


def read_ppm_header(path):
    with open(path, "rb") as f:
        magic = f.readline().strip()
        dims = f.readline().split()
        maxv = f.readline().strip()
    if magic != b"P6" or len(dims) != 2 or maxv != b"255":
        raise ValueError(f"{path}: not a binary P6 PPM")
    return int(dims[0]), int(dims[1])


def load_conf(path):
    conf = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            conf[k.strip()] = v.strip()
    return conf


class Report:
    def __init__(self):
        self.failures = []
        self.lines = []

    def fail(self, msg):
        self.failures.append(msg)
        self.lines.append(f"FAIL: {msg}")

    def ok(self, msg):
        self.lines.append(f"PASS: {msg}")


def check_tier(artdir, tier, conf, rep):
    tdir = os.path.join(artdir, tier)
    metrics_path = os.path.join(tdir, "metrics.json")
    rep.lines.append(f"\n== tier {tier} ==")

    # ---- completeness ---------------------------------------------------
    if not os.path.isfile(metrics_path):
        rep.fail(f"{tier}: missing {metrics_path} (tier not run)")
        return None
    try:
        m = json.load(open(metrics_path))
    except (OSError, ValueError) as e:
        rep.fail(f"{tier}: unreadable metrics.json: {e}")
        return None

    required = ["tier", "params", "budget", "paths", "total",
                "peak_rss_kib", "all_ok", "combined_image_fnv1a64"]
    missing = [k for k in required if k not in m]
    if missing:
        rep.fail(f"{tier}: metrics.json missing keys {missing}")
        return None
    if m["tier"] != tier:
        rep.fail(f"{tier}: metrics reports tier {m['tier']!r}")

    p = m["params"]
    w, h = int(p["width"]), int(p["height"])
    derived = {
        "bvh": w * h * int(p["bvh_spp"]),
        "sample": w * h * int(p["sample_draws"]),
        "material": w * h * int(p["material_spp"]),
    }

    ppm_hashes = {}
    for name in PATHS:
        ppm = os.path.join(tdir, f"{name}.ppm")
        if not os.path.isfile(ppm):
            rep.fail(f"{tier}: missing dataset {name}.ppm")
            continue
        pw, ph = read_ppm_header(ppm)
        if (pw, ph) != (w, h):
            rep.fail(f"{tier}: {name}.ppm is {pw}x{ph}, params say {w}x{h}")
        ppm_hashes[name] = fnv1a64(ppm)
    if len(ppm_hashes) != len(PATHS):
        rep.fail(f"{tier}: incomplete image set, cannot certify tier")
        return None

    # ---- budget ----------------------------------------------------------
    b = m["budget"]
    budget_conf_keys = {
        "bvh": "max_work_units_bvh",
        "sample": "max_work_units_sample",
        "material": "max_work_units_material",
    }
    for name in PATHS:
        entry = m["paths"].get(name)
        if entry is None:
            rep.fail(f"{tier}: metrics has no entry for {name}")
            continue

        conf_bound = int(conf[f"{tier}.{budget_conf_keys[name]}"])
        json_bound = int(b[budget_conf_keys[name]])
        if conf_bound != json_bound:
            rep.fail(f"{tier}/{name}: budget.conf bound {conf_bound} != "
                     f"metrics bound {json_bound}")
        if derived[name] > conf_bound:
            rep.fail(f"{tier}/{name}: scheduled work {derived[name]} exceeds "
                     f"bound {conf_bound}")
        elif derived[name] != conf_bound:
            # The work budget is a scheduled-work budget, so it must equal
            # width*height*spp exactly.  A looser bound would let a path do
            # extra work and still pass.
            rep.fail(f"{tier}/{name}: work bound {conf_bound} != "
                     f"w*h*spp {derived[name]}")
        if int(entry["scheduled_work"]) != derived[name]:
            rep.fail(f"{tier}/{name}: scheduled_work "
                     f"{entry['scheduled_work']} != w*h*spp {derived[name]}")
        if int(entry["measured_work"]) != derived[name]:
            rep.fail(f"{tier}/{name}: measured_work "
                     f"{entry['measured_work']} != scheduled {derived[name]}")
        rep.ok(f"{tier}/{name}: work {entry['measured_work']}/{conf_bound} "
               f"units (scheduled {derived[name]})")

        max_path_ms = float(conf[f"{tier}.max_path_ms"])
        if float(entry["ms_median"]) > max_path_ms:
            rep.fail(f"{tier}/{name}: median {entry['ms_median']:.3f} ms > "
                     f"bound {max_path_ms:.3f} ms")
        rep.ok(f"{tier}/{name}: median {entry['ms_median']:.3f} ms <= "
               f"{max_path_ms:.3f} ms")

        if entry["image_fnv1a64"] != f"{ppm_hashes[name]:016x}":
            rep.fail(f"{tier}/{name}: recorded hash {entry['image_fnv1a64']} "
                     f"!= recomputed {ppm_hashes[name]:016x}")
        rep.ok(f"{tier}/{name}: image hash {ppm_hashes[name]:016x} verified")

    max_total_ms = float(conf[f"{tier}.max_total_ms"])
    total_med = float(m["total"]["ms_median"])
    if total_med > max_total_ms:
        rep.fail(f"{tier}: total median {total_med:.3f} ms > bound "
                 f"{max_total_ms:.3f} ms")
    rep.ok(f"{tier}: total median {total_med:.3f} ms <= {max_total_ms:.3f} ms")

    max_rss = int(conf[f"{tier}.max_peak_rss_kib"])
    rss = int(m["peak_rss_kib"])
    if rss > max_rss:
        rep.fail(f"{tier}: peak RSS {rss} KiB > bound {max_rss} KiB")
    rep.ok(f"{tier}: peak RSS {rss} KiB <= {max_rss} KiB")

    if not m["all_ok"]:
        rep.fail(f"{tier}: harness self-reported all_ok=false")

    return {"metrics": m, "ppm_hashes": ppm_hashes, "w": w, "h": h}


def check_nocheck(nocheck_map, results, rep):
    """PPMs rendered with --no-budget-check must be byte-identical."""
    if not nocheck_map:
        return
    for tier, ndir in nocheck_map.items():
        if tier not in results:
            rep.fail(f"{tier}: --no-budget-check given but checked tier missing")
            continue
        for name in PATHS:
            a = os.path.join(ndir, f"{name}.ppm")
            if not os.path.isfile(a):
                rep.fail(f"{tier}/{name}: missing no-check PPM {a}")
                continue
            ndir_hash = fnv1a64(a)
            ref = results[tier]["ppm_hashes"][name]
            if ndir_hash != ref:
                rep.fail(f"{tier}/{name}: no-check hash {ndir_hash:016x} != "
                         f"checked {ref:016x} (pixels changed when budget "
                         f"checks were disabled)")
            else:
                rep.ok(f"{tier}/{name}: no-check pixels identical "
                       f"({ref:016x})")


def write_report(artdir, results, rep):
    out = os.path.join(artdir, "REPORT.md")
    L = []
    L.append("# Performance budget acceptance -- two tiers\n")
    L.append("Generated by `perf/tools/verify_tiers.py`. Every number is a "
             "measurement from the run below; every bound comes from "
             "`perf/budget.conf`.\n")
    L.append("## Timing table\n")
    L.append("| tier | path | scheduled work | measured work | min ms | "
             "median ms | max ms | median bound ms |")
    L.append("|---|---|---:|---:|---:|---:|---:|---:|")
    for tier in TIERS:
        r = results.get(tier)
        if not r:
            continue
        m = r["metrics"]
        for name in PATHS:
            e = m["paths"][name]
            L.append(f"| {tier} | {name} | {e['scheduled_work']} | "
                     f"{e['measured_work']} | {e['ms_min']:.3f} | "
                     f"{e['ms_median']:.3f} | {e['ms_max']:.3f} | "
                     f"{m['budget']['max_path_ms']:.3f} |")
        t = m["total"]
        L.append(f"| {tier} | **total** | | | {t['ms_min']:.3f} | "
                 f"{t['ms_median']:.3f} | {t['ms_max']:.3f} | "
                 f"{m['budget']['max_total_ms']:.3f} |")
    L.append("")
    L.append("## Peak memory\n")
    L.append("| tier | peak RSS KiB | bound KiB |")
    L.append("|---|---:|---:|")
    for tier in TIERS:
        r = results.get(tier)
        if not r:
            continue
        m = r["metrics"]
        L.append(f"| {tier} | {m['peak_rss_kib']} | "
                 f"{m['budget']['max_peak_rss_kib']} |")
    L.append("")
    L.append("## Image hashes (FNV-1a64 over exact PPM bytes)\n")
    L.append("| tier | bvh | sample | material | combined |")
    L.append("|---|---|---|---|---|")
    for tier in TIERS:
        r = results.get(tier)
        if not r:
            continue
        h = r["ppm_hashes"]
        L.append(f"| {tier} | `{h['bvh']:016x}` | `{h['sample']:016x}` | "
                 f"`{h['material']:016x}` | "
                 f"`{int(r['metrics']['combined_image_fnv1a64'], 16):016x}` |")
    L.append("")
    L.append("## Reproduce\n")
    L.append("```sh")
    L.append("make -C perf accept        # selftest + both tiers + verify + off-switch")
    L.append("")
    L.append("# step by step:")
    L.append("make -C perf clean test")
    L.append("make -C perf low high          # writes artifacts/perf/{low,high}/")
    L.append("python3 -I perf/tools/verify_tiers.py --dir artifacts/perf")
    L.append("# off-switch invariance:")
    L.append("make -C perf no-check")
    L.append("python3 -I perf/tools/verify_tiers.py --dir artifacts/perf \\")
    L.append("        --nocheck-low artifacts/perf/nocheck-low \\")
    L.append("        --nocheck-high artifacts/perf/nocheck-high")
    L.append("```")
    L.append("")
    try:
        with open(out, "w") as f:
            f.write("\n".join(L) + "\n")
        rep.lines.append(f"\nwrote {out}")
    except OSError as e:
        rep.fail(f"cannot write {out}: {e}")


def main(argv):
    artdir = "artifacts/perf"
    budget = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "budget.conf")
    nocheck = {}
    i = 1
    while i < len(argv):
        a = argv[i]
        if a == "--dir" and i + 1 < len(argv):
            artdir = argv[i + 1]
            i += 2
        elif a == "--budget" and i + 1 < len(argv):
            budget = argv[i + 1]
            i += 2
        elif a == "--nocheck-low" and i + 1 < len(argv):
            nocheck["low"] = argv[i + 1]
            i += 2
        elif a == "--nocheck-high" and i + 1 < len(argv):
            nocheck["high"] = argv[i + 1]
            i += 2
        else:
            print(f"usage: {argv[0]} [--dir ART] [--budget FILE] "
                  f"[--nocheck-low DIR] [--nocheck-high DIR]", file=sys.stderr)
            return 2

    if not os.path.isfile(budget):
        print(f"FAIL: budget file not found: {budget}", file=sys.stderr)
        return 2
    conf = load_conf(budget)
    rep = Report()
    results = {}
    for tier in TIERS:
        r = check_tier(os.path.abspath(artdir), tier, conf, rep)
        if r is not None:
            results[tier] = r
        else:
            rep.fail(f"{tier}: incomplete -- tier is not acceptable")

    check_nocheck(nocheck, results, rep)

    # Both tiers must be present before the report can claim acceptance.
    if all(t in results for t in TIERS):
        write_report(os.path.abspath(artdir), results, rep)

    print("\n".join(rep.lines))
    print("\n== summary ==")
    if len(results) != len(TIERS):
        print(f"missing tiers: {[t for t in TIERS if t not in results]}")
    if rep.failures:
        for msg in rep.failures:
            print(f"FAIL: {msg}")
        print("RESULT: FAIL")
        return 1
    print("both tiers complete, every artifact hash verified, "
          "all work/time/memory bounds respected")
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
