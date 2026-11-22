"""
make_figs.py — 生成分析图 (PNG, 150dpi)。
Fig1 MSE 分解与 U 型权衡  Fig2 上界有效性  Fig3 三种最优幅度对比  Fig4 偏差距离。
颜色取 dataviz 参考调色板；中文用 Noto Sans CJK SC。
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

import mse

FONT = "/home/ys/.local/lib/python3.12/site-packages/mplfonts/fonts/NotoSansCJKsc-Regular.otf"
font_manager.fontManager.addfont(FONT)
CJK = font_manager.FontProperties(fname=FONT).get_name()
plt.rcParams["font.family"] = CJK
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 11
plt.rcParams["axes.edgecolor"] = "#c3c2b7"
plt.rcParams["axes.labelcolor"] = "#0b0b0b"
plt.rcParams["xtick.color"] = "#52514e"
plt.rcParams["ytick.color"] = "#52514e"
plt.rcParams["grid.color"] = "#e1e0d9"
plt.rcParams["grid.linewidth"] = 0.8

BLUE, ORANGE, AQUA, VIOLET, RED = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e34948"
INK2, MUTED = "#52514e", "#898781"

CASES = [
    ("gaussian_sigma1.5", "高斯 σ=1.5（窄带）", mse.gaussian_spectrum(1.5)),
    ("gaussian_sigma2.5", "高斯 σ=2.5（中带宽）", mse.gaussian_spectrum(2.5)),
    ("gaussian_sigma4.0", "高斯 σ=4.0（宽带）", mse.gaussian_spectrum(4.0)),
    ("lorentz_p2_w0=2", "洛伦兹 p=2", mse.lorentz_spectrum(2)),
]


def fig1_decomposition():
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    a = np.linspace(0, 0.5, 400)
    for ax, (key, title, S) in zip(axes.ravel(), CASES):
        N, A = mse.mse_components(a, S)
        ax.plot(a, N, color=ORANGE, lw=2, label="N(a) 带内抖动噪声（闪）")
        ax.plot(a, A, color=BLUE, lw=2, label="A(a) 残余混叠（锯齿）")
        ax.plot(a, N + A, color=INK2, lw=2.2, ls="--", label="MSE=N+A")
        a_exh, _, _, _ = mse.exhaustive(S)
        m_star = np.atleast_1d(mse.mse(np.array([a_exh]), S))[0]
        ax.scatter([a_exh], [m_star], s=55, color=RED, zorder=5,
                   label=f"最优 a*={a_exh:.3f}")
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("抖动幅度 a（格距 T=1）")
        ax.set_ylabel("归一化每样本 MSE")
        ax.grid(True, alpha=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.015))
    fig.suptitle("图1  抖动的噪声–混叠权衡：N 随幅度上升、A 随幅度下降，MSE 呈 U 型",
                 fontsize=13.5, y=0.98)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig("figs/fig1_decomposition.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig2_envelope():
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    a = np.linspace(0.02, 0.5, 400)
    for ax, (key, title, S) in zip(
            axes, [CASES[0], CASES[1], CASES[2]]):
        N, A = mse.mse_components(a, S)
        U_N, U_A, mu = mse.envelope_bounds(a, S)
        ax.plot(a, N, color=ORANGE, lw=2, label="N(a) 实际")
        ax.plot(a, U_N, color=ORANGE, lw=1.6, ls=":", label="上界 m2_B·a²/3")
        ax.plot(a, A, color=BLUE, lw=2, label="A(a) 实际")
        ax.plot(a, U_A, color=BLUE, lw=1.6, ls=":", label="上界 A0·min(1,1/(πa)²)")
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("抖动幅度 a")
        ax.set_ylabel("MSE 分量")
        ax.grid(True, alpha=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("图2  采样核频谱推出的误差上界（虚线）处处不低于实际分量（实线）",
                 fontsize=13.5, y=1.0)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    fig.savefig("figs/fig2_envelope.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig3_optima():
    summary = json.load(open("results/summary.json"))
    names = ["gaussian_sigma1.5", "gaussian_sigma2.0", "gaussian_sigma2.5",
             "gaussian_sigma4.0", "lorentz_p2_w0=2", "lorentz_p3_w0=2"]
    labels = ["高斯\nσ1.5", "高斯\nσ2.0", "高斯\nσ2.5", "高斯\nσ4.0",
              "洛伦兹\np2", "洛伦兹\np3"]
    x = np.arange(len(names))
    w = 0.2
    exh = [summary[n]["exhaustive"]["a_star"] for n in names]
    root = [summary[n]["exact_root"]["a_star"] for n in names]
    env = [min(summary[n]["envelope"]["a_star"], 0.5) for n in names]
    qua = [summary[n]["quartic"]["a_star"] or 0.0 for n in names]
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(x - 1.5 * w, exh, w, color=INK2, label="穷举全局最优 a*")
    ax.bar(x - 0.5 * w, root, w, color=RED, label="解析驻点 M′(a)=0")
    ax.bar(x + 0.5 * w, env, w, color=BLUE, label="上界最小化 a_U*（截到[0,½]）")
    ax.bar(x + 1.5 * w, qua, w, color=AQUA, label="四阶闭式 a4*")
    ax.axhline(0.5, color=MUTED, lw=1, ls="--")
    ax.text(len(names) - 0.5, 0.505, "满抖动边界 a=½", color=MUTED, fontsize=9,
            ha="right", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("最优抖动幅度 a")
    ax.set_ylim(0, 0.62)
    ax.set_title("图3  三种解析最优幅度与穷举最优的对比（解析驻点与穷举在所有用例重合）")
    ax.grid(True, axis="y", alpha=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig("figs/fig3_optima.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig4_distance():
    summary = json.load(open("results/summary.json"))
    names = ["gaussian_sigma2.5", "gaussian_sigma4.0", "lorentz_p2_w0=2"]
    labels = ["高斯 σ2.5", "高斯 σ4.0", "洛伦兹 p2"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))

    # 左：解析驻点与穷举的实测距离（对数轴，含网格步参考线）
    d_root = [summary[n]["exact_root"]["d_measured"] for n in names]
    h = summary[names[0]]["exhaustive"]["grid_step"]
    ax = axes[0]
    ax.bar(labels, d_root, color=RED, width=0.55, label="实测 |â−a*|")
    ax.axhline(h, color=MUTED, ls="--", lw=1.2, label=f"穷举网格步 h={h:.5f}")
    for i, v in enumerate(d_root):
        ax.text(i, v + 2e-5, f"{v:.1e}", ha="center", fontsize=10, color=INK2)
    ax.set_ylabel("距离 |â − a*|")
    ax.set_title("解析驻点 vs 穷举：距离 ≤ 1 个网格步")
    ax.set_ylim(0, max(h * 1.6, max(d_root) * 1.4))
    ax.grid(True, axis="y", alpha=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9.5)

    # 右：两个闭式近似的实测距离与强凸上界
    ax = axes[1]
    x = np.arange(len(names))
    wbar = 0.35
    d_env = [summary[n]["envelope"]["d_measured"] for n in names]
    u_env = [summary[n]["envelope"]["d_upper_convex"] for n in names]
    d_qua = [summary[n]["quartic"]["d_measured"] for n in names]
    u_qua = [summary[n]["quartic"]["d_upper_convex"] for n in names]
    ax.bar(x - wbar / 2, d_env, wbar, color=BLUE, label="上界法 实测距离")
    ax.bar(x + wbar / 2, u_env, wbar, color=BLUE, alpha=0.4, label="上界法 强凸距离上界")
    ax.bar(x - wbar / 2, d_qua, wbar * 0.5, color=AQUA, label="四阶法 实测距离")
    ax.bar(x + wbar / 2, u_qua, wbar * 0.5, color=AQUA, alpha=0.4,
           label="四阶法 强凸距离上界")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("距离 / 距离上界")
    ax.set_title("闭式近似的偏差与强凸上界 d≤√(2Δ/κ)+h")
    ax.grid(True, axis="y", alpha=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8.5, ncol=2)
    fig.suptitle("图4  偏差上界核验：实测距离始终落在强凸上界之内", fontsize=13.5, y=1.02)
    fig.tight_layout()
    fig.savefig("figs/fig4_distance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig1_decomposition()
    fig2_envelope()
    fig3_optima()
    fig4_distance()
    print("wrote figs/fig1..fig4.png")
