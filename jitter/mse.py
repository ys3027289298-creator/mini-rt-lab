"""
mse.py — 均匀抖动(分层)采样的重构 MSE 模型、误差上界与最优抖动幅度。

约定
----
* 采样间距 T = 1（频率单位：弧度/采样），基带 B = [-π, π]，奈奎斯特频率 π。
* 零均值平稳信号 x(t)，双边功率谱 S(ω)，归一化 E = (1/2π)∫S dω = 1。
* 抖动样本 t_n = n + U_n，U_n iid ~ Uniform(-a, a)，0 ≤ a ≤ 1/2。
  a = 0 退化为确定性格点采样；a = 1/2 为满抖动(每格内全覆盖)。
* 抖动特征函数(采样核频谱) Φ_U(ω) = E[e^{jωU}] = sinc(aω)，sinc(x)=sin(x)/x。

Balakrishnan (1962) 框架下，理想低通重构的每样本基带 MSE 分解为
  MSE(a) = N(a) + A(a)
  N(a) = 1/(2π) ∫_B          [1 - sinc²(aω)] S(ω) dω              带内抖动噪声
  A(a) = 1/(2π) ∫_B Σ_{k≠0}       sinc²(a(ω+2πk)) S(ω+2πk) dω    混叠(相干→被核压制)
二者对 a 单调性相反：N 随 a 增大(画面"闪")，A 随 a 减小(锯齿/摩尔纹消失)。
"""

import json
import numpy as np

# ---------------------------------------------------------------- 频率网格
# 令 dw = 2π/n_per，使混叠平移 2π 恰好对应整数个网格点，混叠求和退化为数组平移。
N_PER = 4096          # 每个 2π 周期内的网格点数
N_PERIODS = 16        # 网格两侧各覆盖多少个 2π 周期
HALF = N_PER // 2     # 基带 B=[-π,π] 对应的半宽索引


def freq_grid():
    n = N_PER * N_PERIODS
    j = np.arange(-n, n + 1)
    w = j * (2.0 * np.pi / N_PER)
    band = np.abs(j) <= HALF
    return w, j, band


def sinc(x):
    """未归一化 sinc: sin(x)/x，x=0 处取 1。"""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    nz = np.abs(x) > 1e-14
    out[nz] = np.sin(x[nz]) / x[nz]
    return out


def d_sinc_da(a, w):
    """d/da sinc(a w)，a→0 时用 -a w²/3 极限，全部以数组返回。"""
    x = a * w
    out = np.empty_like(x)
    small = np.abs(x) < 1e-6
    out[small] = -a * w[small] ** 2 / 3.0
    nx = x[~small]
    out[~small] = w[~small] * (nx * np.cos(nx) - np.sin(nx)) / nx ** 2
    return out


# ---------------------------------------------------------------- 信号频谱
def gaussian_spectrum(sigma):
    """高斯谱 S(ω) ∝ exp(-ω²/(2σ²))，闭合式二阶矩 m2 = σ²。"""
    w, _, _ = freq_grid()
    S = np.exp(-(w ** 2) / (2.0 * sigma ** 2))
    S /= np.trapezoid(S, w) / (2.0 * np.pi)
    return S


def lorentz_spectrum(p, w0=2.0):
    """广义洛伦兹谱 S(ω) ∝ 1/(1+(ω/w0)²)^p；m2 = w0²/(2p-3) (p>3/2)。"""
    w, _, _ = freq_grid()
    S = 1.0 / (1.0 + (w / w0) ** 2) ** p
    S /= np.trapezoid(S, w) / (2.0 * np.pi)
    return S


# ---------------------------------------------------------------- 混叠索引
# 频率网格的位置下标：j = arange(-n,n+1)，位置 p = j + n。
# 对基带内每个位置 p_b，收集它的全部混叠位置 p_b + k*N_PER (k≠0)。
def alias_index_matrix():
    _, j, band = freq_grid()
    n = N_PER * N_PERIODS
    n_total = j.size
    pos = j + n                                   # 位置下标 0..2n
    pb = pos[band]                                # 基带位置
    ks = np.arange(-N_PERIODS, N_PERIODS + 1)
    ks = ks[ks != 0]
    idx = pb[:, None] + ks[None, :] * N_PER
    valid = (idx >= 0) & (idx < n_total)
    idx = np.where(valid, idx, 0)                 # 越界位用 0 占位，乘 0 屏蔽
    return idx, valid, pb


_ALIAS_CACHE = None


def _alias():
    global _ALIAS_CACHE
    if _ALIAS_CACHE is None:
        _ALIAS_CACHE = alias_index_matrix()
    return _ALIAS_CACHE


# ---------------------------------------------------------------- MSE 分量
def mse_components(a, S):
    """返回 (N(a), A(a))，支持 a 为标量或数组。"""
    w, _, band = freq_grid()
    idx, valid, pb = _alias()
    a = np.atleast_1d(np.asarray(a, dtype=float))
    N_out = np.empty(a.size)
    A_out = np.empty(a.size)
    wb = w[band]
    Sa = S[idx] * valid                      # 混叠位置上的谱
    wa = w[idx]                              # 混叠频率 u = ω+2πk
    chunk = 64
    for s in range(0, a.size, chunk):
        aa = a[s:s + chunk].reshape(-1, 1)                 # (c,1)
        q_b = sinc(aa * wb[None, :]) ** 2                  # (c, nband)
        N_out[s:s + chunk] = np.trapezoid(
            (1.0 - q_b) * S[None, pb], wb, axis=1) / (2.0 * np.pi)
        q_a = sinc(aa[:, :, None] * wa[None, :, :]) ** 2   # (c, nband, nk)
        integrand = (q_a * Sa[None, :, :]).sum(axis=2)
        A_out[s:s + chunk] = np.trapezoid(integrand, wb, axis=1) / (2.0 * np.pi)
    return N_out, A_out


def mse(a, S):
    n, al = mse_components(a, S)
    return n + al


def dmse_da(a, S):
    """MSE 的解析导数 d/da [N+A]，用于驻点方程 M'(a)=0。"""
    w, _, band = freq_grid()
    idx, valid, pb = _alias()
    wb = w[band]
    Sa = S[idx] * valid
    wa = w[idx]
    a = float(a)
    # d/da sinc²(aω) = 2 sinc(aω) * d/da sinc(aω)；向量化一次算全部频率
    xb = a * wb
    dq_b = 2.0 * sinc(xb) * d_sinc_da(a, wb)
    dN = np.trapezoid(-dq_b * S[pb], wb) / (2.0 * np.pi)
    xa = a * wa
    dq_a = np.where(valid, 2.0 * sinc(xa) * d_sinc_da(a, wa), 0.0)
    dA = np.trapezoid((dq_a * Sa).sum(axis=1), wb) / (2.0 * np.pi)
    return dN + dA


# ---------------------------------------------------------------- 谱矩与上界
def spectral_moments(S):
    """计算闭合式最优幅度所需的各阶谱矩(全部按 1/2π 归一化)。"""
    w, _, band = freq_grid()
    idx, valid, _ = _alias()
    E = np.trapezoid(S, w) / (2.0 * np.pi)
    m2 = np.trapezoid(w ** 2 * S, w) / (2.0 * np.pi)
    m4 = np.trapezoid(w ** 4 * S, w) / (2.0 * np.pi)
    m2_B = np.trapezoid(w[band] ** 2 * S[band], w[band]) / (2.0 * np.pi)
    m4_B = np.trapezoid(w[band] ** 4 * S[band], w[band]) / (2.0 * np.pi)
    # 确定性采样(a=0)的混叠能量 A0
    A0 = np.trapezoid((S[idx] * valid).sum(axis=1), w[band]) / (2.0 * np.pi)
    return dict(E=E, m2=m2, m4=m4, m2_B=m2_B, m4_B=m4_B, A0=A0)


def envelope_bounds(a, S):
    """
    由采样核频谱给出的两条上界(对全部 a≥0 成立):
      N(a) = ∫(1-sinc²)S ≤ a²/3 · m2_B          （1-sinc²(x) ≤ x²/3）
      A(a) ≤ A0 · min(1, 1/(πa)²)                （混叠频率 |u|≥π，sinc²(x)≤min(1,1/x²)）
    总上界 U(a) = m2_B a²/3 + A0/(π²a²)（在 a≥1/π 的 1/a² 分支上）。
    """
    mu = spectral_moments(S)
    a = np.atleast_1d(np.asarray(a, dtype=float))
    U_N = mu["m2_B"] * a ** 2 / 3.0
    U_A = mu["A0"] * np.minimum(1.0, 1.0 / (np.pi * a) ** 2)
    return U_N, U_A, mu


def envelope_optimum(S):
    """
    最小化总上界 U(a)=m2_B a²/3 + A0/(π²a²)：
      a_U* = (3 A0 / (π² m2_B))^{1/4}
    返回 (a_star, 状态)。状态说明分支条件：
      * 若 a_U* < 1/π，1/a² 分支尚未生效，上界在饱和区单调增 → 上界法给不出内点，
        只能判定 a*=0 一侧，需用精确驻点方程；
      * 若 a_U* > 1/2，可行域内 U 单调减 → 上界推荐满抖动边界；
      * 否则内点有效，且 U''=2m2_B/3 + 6A0/(π²a⁴)>0 保证唯一最小。
    """
    mu = spectral_moments(S)
    a_u = (3.0 * mu["A0"] / (np.pi ** 2 * mu["m2_B"])) ** 0.25
    if a_u < 1.0 / np.pi:
        status = "sat:envelope_too_loose_use_exact_root"
    elif a_u > 0.5:
        status = "boundary:full_jitter"
    else:
        status = "interior"
    return a_u, status, mu


def quartic_optimum(S):
    """
    小抖动四阶展开 MSE(a) ≈ c2 a² + c4 a⁴ 的闭式最优幅度：
      c2 = (2 m2_B - m2)/3，  c4 = (m4 - 2 m4_B)/90
      a4* = sqrt(-c2/(2 c4)) = sqrt(30 (m2 - 2m2_B)/(m4 - 2m4_B))
    条件: c2<0（带外二阶矩占主导，否则 a*=0）且 c4>0（四阶项把曲线拉回）。
    """
    mu = spectral_moments(S)
    c2 = (2.0 * mu["m2_B"] - mu["m2"]) / 3.0
    c4 = (mu["m4"] - 2.0 * mu["m4_B"]) / 90.0
    if c2 >= 0:
        return 0.0, "c2>=0:no_interior_minimum_a*=0", c2, c4, mu
    if c4 <= 0:
        return float("nan"), "c4<=0:quartic_unbounded_need_higher_order", c2, c4, mu
    return float(np.sqrt(-c2 / (2.0 * c4))), "interior", c2, c4, mu


# ---------------------------------------------------------------- 精确驻点与穷举
def exact_root(S):
    """
    全局最优：MSE 可以是双峰曲线(内点驻点与边界 a=1/2 几乎简并)，
    因此枚举驻点方程 M'(a)=0 的*全部*根，连同两端点 a=0、a=1/2 一起比较 MSE。
    每个符号翻转换 Brent 精化；并报告 M''(a*) 符号(局部最小条件)。
    """
    from scipy import optimize
    grid = np.linspace(0.0, 0.5, 4001)
    g = np.array([dmse_da(aa, S) for aa in grid])
    candidates = [0.0, 0.5]
    for i in np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]:
        candidates.append(
            optimize.brentq(dmse_da, grid[i], grid[i + 1], args=(S,), xtol=1e-12))
    vals = [float(np.atleast_1d(mse(a, S))[0]) for a in candidates]
    k = int(np.argmin(vals))
    a_star = float(candidates[k])
    if a_star in (0.0, 0.5):
        return a_star, f"boundary:a={a_star};n_stationary={len(candidates)-2}"
    h = 1e-4
    d2 = (dmse_da(a_star + h, S) - dmse_da(a_star - h, S)) / (2 * h)
    status = (f"interior_convex:n_stationary={len(candidates)-2}:M2={d2:.3e}"
              if d2 > 0 else f"interior_nonconvex:M2={d2:.3e}")
    return a_star, status


def exhaustive(S, n_grid=4001):
    """在 a∈[0,1/2] 细网格上穷举，返回 (最优点, 网格步长, 全部a, 全部MSE)。"""
    a = np.linspace(0.0, 0.5, n_grid)
    M = mse(a, S)
    i = int(np.argmin(M))
    return float(a[i]), float(a[1] - a[0]), a, M


def gap_distance_bound(a_hat, a_star, S, h_grid):
    """
    解析最优点 â 相对穷举最优点 a* 的偏差界。
      Δ = M(â)-M(a*) ≥ 0。
      * 强凸上界（主结果）：a* 邻域 M''≥κ>0 时 Δ ≥ κ d²/2
        → d ≤ sqrt(2Δ/κ)。κ 取 [min(â,a*), max] 区间上 M'' 的最小值，保证覆盖。
      * Lipschitz 下界（仅核验）：Δ ≤ L·d → d ≥ Δ/L。
      * 穷举自身有 ±h_grid/2 的网格不确定度，上界中计入。
    另外检查 â 是否落在 a* 的一个网格步长内（"落在穷举最优点附近"的直接证据）。
    """
    M_hat = float(np.atleast_1d(mse(a_hat, S))[0])
    M_star = float(np.atleast_1d(mse(a_star, S))[0])
    delta = max(M_hat - M_star, 0.0)
    h = 1e-4
    kappa_star = (dmse_da(a_star + h, S) - dmse_da(a_star - h, S)) / (2 * h)
    lo, hi = sorted((a_hat, a_star))
    probe = np.linspace(max(lo, 1e-3), hi, 30) if hi > 1e-3 else np.array([h])
    kappas = np.array([
        (dmse_da(aa + h, S) - dmse_da(aa - h, S)) / (2 * h) for aa in probe])
    kappa_min = float(np.min(kappas))
    L = float(max(abs(dmse_da(aa, S)) for aa in
                  np.linspace(lo, hi + 1e-12, 50)))
    d_meas = abs(a_hat - a_star)
    d_upper = (float(np.sqrt(2 * delta / kappa_min)) + h_grid
               if (kappa_min > 0 and delta > 0) else float("inf"))
    d_lower = float(delta / L) if L > 0 else 0.0
    return dict(delta=delta, kappa_at_star=float(kappa_star),
                kappa_min_interval=kappa_min, L=L,
                d_measured=d_meas,
                d_upper_convex=d_upper,
                d_lower_lip=d_lower,
                within_one_gridstep=bool(d_meas <= h_grid + 1e-12),
                grid_uncertainty=h_grid)


def run_case(name, S):
    """对一个谱跑完整流程，汇总为可 JSON 序列化的字典。"""
    a_exh, h_grid, a_grid, M_grid = exhaustive(S)
    a_root, root_status = exact_root(S)
    a_u, env_status, mu = envelope_optimum(S)
    a4, q_status, c2, c4, _ = quartic_optimum(S)
    N0, A0v = mse_components(np.array([0.0]), S)
    N_full, A_full = mse_components(np.array([0.5]), S)
    out = dict(
        name=name,
        moments={k: float(v) for k, v in mu.items()},
        c2=float(c2), c4=float(c4),
        mse_a0=float((N0 + A0v)[0]),
        mse_afull=float((N_full + A_full)[0]),
        exhaustive=dict(a_star=a_exh, grid_step=h_grid,
                        mse_star=float(np.min(M_grid))),
        exact_root=dict(a_star=a_root, status=root_status,
                        **gap_distance_bound(a_root, a_exh, S, h_grid)),
        envelope=dict(a_star=float(a_u), status=env_status,
                      **gap_distance_bound(float(np.clip(a_u, 0, 0.5)),
                                           a_exh, S, h_grid)),
        quartic=dict(a_star=float(a4) if np.isfinite(a4) else None,
                     status=q_status,
                     **(gap_distance_bound(a4, a_exh, S, h_grid)
                        if np.isfinite(a4) else dict())),
    )
    return out, a_grid, M_grid


def main():
    cases = {
        "gaussian_sigma1.5": gaussian_spectrum(1.5),
        "gaussian_sigma2.0": gaussian_spectrum(2.0),
        "gaussian_sigma2.5": gaussian_spectrum(2.5),
        "gaussian_sigma4.0": gaussian_spectrum(4.0),
        "lorentz_p2_w0=2": lorentz_spectrum(2),
        "lorentz_p3_w0=2": lorentz_spectrum(3),
    }
    summary = {}
    curves = {}
    for name, S in cases.items():
        info, a_grid, M_grid = run_case(name, S)
        summary[name] = info
        curves[name] = dict(a=a_grid.tolist(), mse=M_grid.tolist())
        print(f"== {name}: exhaust a*={info['exhaustive']['a_star']:.5f} "
              f"root a*={info['exact_root']['a_star']:.5f} "
              f"env a*={info['envelope']['a_star']:.5f} "
              f"quartic a*={info['quartic']['a_star']}")
    with open("results/summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    with open("results/curves.json", "w") as f:
        json.dump(curves, f)
    print("wrote results/summary.json, results/curves.json")


if __name__ == "__main__":
    main()
