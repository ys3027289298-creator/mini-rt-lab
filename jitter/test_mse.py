"""
test_mse.py — 模型与数值实现的独立健全性测试。全部通过才允许出图/写报告。
覆盖:
  1. 归一化与闭合谱矩 (高斯 m2=σ², m4=3σ⁴；洛伦兹 m2=w0²/(2p-3))
  2. a=0 退化: N(0)=0, MSE(0)=A0=确定性采样混叠能量
  3. A0 与独立的谱尾概率一致
  4. N 随 a 单调增(闪)、A 随 a 单调减(锯齿)、均非负
  5. 解析导数 vs 中心有限差分
  6. 解析驻点 vs 穷举全局最优一致到一个网格步
  7. 包络上界确实是上界 (U_N≥N, U_A≥A 处处成立)
  8. 窄带信号最优 a*=0 (确定性采样退化)
"""
import numpy as np
from scipy.stats import norm

import mse

PASS, FAIL = 0, 0


def check(name, cond, detail=""):
    global PASS, FAIL
    ok = bool(cond)
    PASS += ok
    FAIL += (not ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")


def test_moments():
    for sig in [1.5, 2.5, 4.0]:
        S = mse.gaussian_spectrum(sig)
        mu = mse.spectral_moments(S)
        check(f"gauss σ={sig} m2=σ²", abs(mu["m2"] - sig ** 2) < 1e-3,
              f"{mu['m2']:.4f} vs {sig**2:.4f}")
        check(f"gauss σ={sig} m4=3σ⁴", abs(mu["m4"] - 3 * sig ** 4) < 1e-2,
              f"{mu['m4']:.4f} vs {3*sig**4:.4f}")
        check(f"gauss σ={sig} E=1", abs(mu["E"] - 1) < 1e-9)
    p, w0 = 3, 2.0
    mu = mse.spectral_moments(mse.lorentz_spectrum(p, w0))
    check("lorentz p3 m2=w0²/(2p-3)",
          abs(mu["m2"] - w0 ** 2 / (2 * p - 3)) < 1e-2,
          f"{mu['m2']:.4f} vs {w0**2/(2*p-3):.4f}")


def test_zero_jitter_degenerate():
    for sig in [1.5, 2.5, 4.0]:
        S = mse.gaussian_spectrum(sig)
        N, A = mse.mse_components(np.array([0.0]), S)
        mu = mse.spectral_moments(S)
        check(f"σ={sig} N(0)=0", N[0] < 1e-12)
        check(f"σ={sig} MSE(0)=A0", abs(N[0] + A[0] - mu["A0"]) < 1e-12)
        # 独立理论: 高斯谱 |ω|>π 的双边尾概率
        A0_th = 2 * norm.sf(np.pi, scale=sig)
        check(f"σ={sig} A0=尾概率", abs(mu["A0"] - A0_th) < 1e-3,
              f"{mu['A0']:.5f} vs {A0_th:.5f}")


def test_monotonicity():
    a = np.linspace(0, 0.5, 200)
    for sig in [1.5, 2.5, 4.0]:
        N, A = mse.mse_components(a, mse.gaussian_spectrum(sig))
        check(f"σ={sig} N 单调增(抖动噪声)", np.all(np.diff(N) >= -1e-12))
        check(f"σ={sig} A 单调减(混叠被压制)", np.all(np.diff(A) <= 1e-12))
        check(f"σ={sig} N,A 非负", np.all(N >= -1e-12) and np.all(A >= -1e-12))


def test_derivative():
    for sig in [1.5, 2.0, 2.5, 4.0]:
        S = mse.gaussian_spectrum(sig)
        for aa in [1e-2, 0.1, 0.3, 0.49]:
            h = 1e-5
            fd = (np.atleast_1d(mse.mse(aa + h, S))[0]
                  - np.atleast_1d(mse.mse(aa - h, S))[0]) / (2 * h)
            an = mse.dmse_da(aa, S)
            check(f"σ={sig} a={aa} 导数一致",
                  abs(an - fd) <= 1e-3 * max(abs(fd), 1e-12),
                  f"analytic={an:.2e} fd={fd:.2e}")


def test_envelope_is_upper_bound():
    a = np.linspace(0.02, 0.5, 100)
    for sig in [1.5, 2.5, 4.0]:
        S = mse.gaussian_spectrum(sig)
        N, A = mse.mse_components(a, S)
        U_N, U_A, _ = mse.envelope_bounds(a, S)
        check(f"σ={sig} U_N≥N 处处", np.all(U_N >= N - 1e-12))
        check(f"σ={sig} U_A≥A 处处", np.all(U_A >= A - 1e-12))


def test_optima_agree():
    cases = {
        "g1.5": mse.gaussian_spectrum(1.5),
        "g2.0": mse.gaussian_spectrum(2.0),
        "g2.5": mse.gaussian_spectrum(2.5),
        "g4.0": mse.gaussian_spectrum(4.0),
        "lor2": mse.lorentz_spectrum(2),
        "lor3": mse.lorentz_spectrum(3),
    }
    for name, S in cases.items():
        a_exh, h, _, _ = mse.exhaustive(S)
        a_root, status = mse.exact_root(S)
        check(f"{name} 驻点∈穷举一个网格步",
              abs(a_root - a_exh) <= h + 1e-9,
              f"root={a_root:.5f} exhaust={a_exh:.5f} h={h:.5f} [{status}]")


def test_narrowband_prefers_zero():
    for sig in [1.5, 2.0]:
        a_exh, _, _, _ = mse.exhaustive(mse.gaussian_spectrum(sig))
        check(f"窄带 σ={sig} 全局最优 a*=0", a_exh < 1e-9)
    # 四阶判据 c2>=0 与之一致
    for sig in [1.5, 2.0]:
        a4, st, *_ = mse.quartic_optimum(mse.gaussian_spectrum(sig))
        check(f"窄带 σ={sig} 四阶判据给 a*=0", a4 == 0.0 and st.startswith("c2"))


if __name__ == "__main__":
    test_moments()
    test_zero_jitter_degenerate()
    test_monotonicity()
    test_derivative()
    test_envelope_is_upper_bound()
    test_optima_agree()
    test_narrowband_prefers_zero()
    print(f"\n{PASS} passed, {FAIL} failed, all_ok={FAIL == 0}")
    raise SystemExit(1 if FAIL else 0)
