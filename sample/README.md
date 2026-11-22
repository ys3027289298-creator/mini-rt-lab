# sample

`make test` 跑自检。

- `src/sample.c` — 分层与低差异采样、方差估计。
- `src/alloc.c` + `src/alloc.h` — 多光源样本数最优分配：
  - 目标 `min sum_i var_i / n_i`，约束 `sum_i cost_i * n_i <= B`、`n_i >= m`；
    连续最优解 `n_i* ∝ sigma_i / sqrt(cost_i)`（水位法求解，含最小样本约束），
    其方差 `theoretical_var` 是所有可行整数分配的下界。
  - 整数分配对连续解向下取整，再按“单位成本边际方差下降”贪心花完剩余预算。
  - `alloc_bruteforce` 穷举小预算最优值；`alloc_deviation_bound` 给出
    `(allocated_var - theoretical_var) / theoretical_var` 的严格上界，
    因此也覆盖与穷举最优的偏差。
  - 方差缺失（`var_known == 0`、负值、NaN、inf）时用保守值：
    调用者指定的值优先，否则取已知方差最大值（全未知时为 1）；
    通过 `fallback_mask` / `fallback_count` / `fallback_var` 记录。
    已知的 0 方差是合法值（如被完全遮挡），不会触发回退。
