# allocate

按像素方差做自适应样本预算分配，替代每像素固定样本数。`make test` 跑自检
（确定性 LCG，同种子跨进程输出逐字节一致）。

## 解决的问题

固定每像素 `n` 个样本时，高方差像素（边缘、焦散、光泽反射）算不够，平滑
区域的样本被浪费。给定总预算 `B`，像素 `i` 的估计量方差为
`s_i^2 / n_i`，图像均方误差

```
MSE = (1/N) * Σ_i s_i^2 / n_i
```

## 预分配：Neyman 分配

在约束 `Σ n_i = B`、`n_i ≥ n_min` 下最小化 MSE，拉格朗日乘子给出
连续最优解 `n_i ∝ s_i`（标准差，不是方差）。实现要点：

- 每像素保底 `n_min ≥ 2`（少于 2 个样本无法估方差）；低于保底的像素
  钉在保底值，权重从剩余像素重新分配（water-filling 迭代）；
- 整数化用最大余数法，预算恰好花完；平局按当前样本数、再按像素下标决定，
  保证结果唯一可复现；
- 方差估计为零的像素只拿保底值。

输出同时给出每像素样本数和分配比例 `f_i = n_i / B`。

## 与均匀分配的误差比对

同预算 `B` 下均匀分配 `n_i = B/N`：

```
RMSE_uniform  = sqrt(Σ s_i^2 / B)
RMSE_adaptive = sqrt((Σ s_i)^2 / (N B))     # Neyman 连续最优
ratio         = RMSE_adaptive / RMSE_uniform
              = (Σ s_i) / sqrt(N · Σ s_i^2) ≤ 1
```

场景内像素噪声越不均匀，ratio 越小（自检场景 8 像素、σ ∈ {0.5…4} 时
ratio ≈ 0.81；理论下界 1/√N）。自检还用 400 次独立蒙特卡洛实验验证经验
ratio 与预测值偏差 < 0.15，且自适应 MSE 显著低于均匀。

## 目标误差与收敛判据（可复现）

达到目标 RMSE `ε` 所需预算直接由上式解出：

```
B* = (Σ s_i)^2 / (N ε^2)，再向上取整并抬到保底预算
```

`converge()` 的停止规则（固定种子 LCG + Box–Muller，无 `rand()`/`time()`，
同配置重跑结果逐字节相同）：

1. 先导阶段：每像素采 `n_min` 个样本，Welford 在线估 `s_i^2`；
2. 每轮用当前方差做 Neyman 分配，只追采增量（样本数只增不减）；
3. 预测 RMSE = `sqrt(mean_i s_i^2 / n_i)`；
4. **连续两轮** RMSE ≤ 目标才判定收敛（双轮确认，抑制先导噪声误判）；
5. 预算按 `max(2B, B*)` 增长，直到收敛或触及 `max_budget`；
6. 上限内达不到目标返回 `ALLOC_E_TARGET` 并给出所需预算；硬上限
   `MAX_ROUNDS=16` 轮。

## 预算为零 / 非法输入

最小可用预算 = `N · n_min`（每个像素至少 2 个样本）。`budget = 0`、
负数或低于该值一律返回 `ALLOC_E_BUDGET`，并通过出参回报最小可用值：

```
zero_budget: rc=1 min_required=32
```

目标 RMSE ≤ 0、像素数越界等参数错误返回 `ALLOC_E_PARAM`。

## 自检输出（节选）

```
--- allocation demo: budget=512 n_min=4 min_usable=32 ---
pixel  sigma  samples  fraction
    0   0.50       17  0.0332
    ...
    7   4.00      137  0.2676
rmse_uniform=0.288111 rmse_adaptive=0.234377 error_ratio=0.813495
converged=1 rounds=3 final_budget=128 final_rmse=0.487229 target=0.800
all_ok=1
```

## 接入真实渲染器

当前 `pixel_sample()` 是合成像素（正态噪声，代替一次相机光线追踪）。
接入时把它替换为逐像素的路径追踪求值，把先导/追采循环挂到渲染主循环，
其余分配与收敛逻辑不变。
