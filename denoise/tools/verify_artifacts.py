#!/usr/bin/env python3
"""独立核验 denoise 产物：从 CSV / f32 / PGM 重算所有结论，不信任 C 自检输出。

用法: python3 -I tools/verify_artifacts.py <artifacts_dir>
"""
import csv
import struct
import sys
import math

W, H, PX = 64, 48, 64 * 48
WARMUP, FLICKER_EPS, FRAMES = 12, 2e-3, 40
MAGIC = 0x52544731


def read_f32(path):
    with open(path, "rb") as f:
        data = f.read()
    magic, w, h = struct.unpack_from("<Iii", data, 0)
    assert magic == MAGIC, f"{path}: bad magic {magic:#x}"
    assert (w, h) == (W, H), f"{path}: bad dims {w}x{h}"
    n = len(data[12:]) // 4
    return list(struct.unpack(f"<{n}f", data[12:]))


def read_pgm(path):
    with open(path, "rb") as f:
        assert f.readline().strip() == b"P5"
        w, h = map(int, f.readline().split())
        maxv = int(f.readline())
        pix = f.read()
    assert (w, h) == (W, H) and maxv == 255 and len(pix) == PX
    return pix


def main():
    art = sys.argv[1]
    fails = []

    def expect(name, ok, detail=""):
        print(f"verify_{name}={'PASS' if ok else 'FAIL'} {detail}".rstrip())
        if not ok:
            fails.append(name)

    rows = list(csv.DictReader(open(f"{art}/variance_stats.csv")))
    assert len(rows) == FRAMES
    lm = [float(r["legacy_mean"]) for r in rows]
    um = [float(r["unified_mean"]) for r in rows]
    umv = [float(r["unified_mean_var"]) for r in rows]
    lmse = [float(r["legacy_mse"]) for r in rows]
    umse = [float(r["unified_mse"]) for r in rows]

    def flicker(means):
        hits, worst = [], 0.0
        for t in range(WARMUP + 1, FRAMES):
            rel = abs(means[t] - means[t - 1]) / max(abs(means[t - 1]), 1e-9)
            worst = max(worst, rel)
            if rel > FLICKER_EPS:
                hits.append((t, rel))
        return len(hits), worst

    lf, lworst = flicker(lm)
    uf, uworst = flicker(um)
    expect("legacy_flicker_present", lf > 0, f"frames={lf} worst_rel={lworst:.6f}")
    expect("unified_flicker_zero", uf == 0, f"frames={uf} worst_rel={uworst:.6f}")

    le = sum(lm[t] for t in range(WARMUP, FRAMES) if t % 2 == 0)
    lo = sum(lm[t] for t in range(WARMUP, FRAMES) if t % 2 == 1)
    ue = sum(um[t] for t in range(WARMUP, FRAMES) if t % 2 == 0)
    uo = sum(um[t] for t in range(WARMUP, FRAMES) if t % 2 == 1)
    ratio = (le / 14) / (lo / 14)
    uratio = (ue / 14) / (uo / 14)
    expect("legacy_even_odd_ratio", ratio > 1.05, f"ratio={ratio:.4f}")
    expect("unified_even_odd_ratio", abs(uratio - 1.0) < 0.002, f"ratio={uratio:.4f}")

    # 旁路：输入与输出 float 载荷必须逐字节一致
    binp = read_f32(f"{art}/bypass_input.f32")
    bout = read_f32(f"{art}/bypass_output.f32")
    expect("bypass_byte_identical",
           len(binp) == PX and binp == bout, f"pixels={len(binp)}")

    # 方差：均值方差应单调下降（允许 1e-12 数值噪声），末帧 < 首帧 1/4
    mono = all(umv[t + 1] <= umv[t] + 1e-12 for t in range(FRAMES - 1))
    expect("mean_var_monotone", mono)
    expect("mean_var_quarter", umv[-1] < 0.25 * umv[0],
           f"first={umv[0]:.6f} final={umv[-1]:.6f} ratio={umv[-1]/umv[0]:.4f}")
    expect("unified_mse_better", umse[-1] < lmse[-1],
           f"legacy={lmse[-1]:.6f} unified={umse[-1]:.6f}")

    # golden：旧路径输出冻结文件存在且形状正确
    golden = read_f32(f"{art}/golden_legacy_off.f32")
    expect("golden_shape", len(golden) == PX)

    # 噪点图像摘要：PGM 可读、统计量合理（有噪但非全黑/全白）
    noisy = read_pgm(f"{art}/noisy_frame00.pgm")
    nmean = sum(noisy) / PX
    nvar = sum((p - nmean) ** 2 for p in noisy) / PX
    expect("noisy_pgm_sane", 5 < nmean < 250 and nvar > 25,
           f"mean={nmean:.2f} std={math.sqrt(nvar):.2f}")
    for name in ("legacy_output_frame39.pgm", "unified_output_frame39.pgm"):
        read_pgm(f"{art}/{name}")  # 仅校验可解析与尺寸
    print(f"noisy_summary pixels={PX} mean={nmean:.2f} std={math.sqrt(nvar):.2f} "
          f"min={min(noisy)} max={max(noisy)}")

    print(f"OVERALL={'PASS' if not fails else 'FAIL'} failed={fails}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
