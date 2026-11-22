# sample

`make test` 跑自检。

除分层/低差异采样与方差估计外，包含光子发射功率补偿：按光源功率 `P_i` 与目标发射分布 `q_i` 推每光子的发射权重 `w_i = (P_i / P_total) / q_i`，估计量 `X = c_i * w_i` 对 `sum_i c_i P_i / P_total` 无偏；与均匀发射（`q_i = 1/n`）比对方差，测试输出 `var_ratio = var_target / var_uniform`。

约束：目标分布必须积分为一（`sum q_i == 1`），权重为负一律拒绝（含负功率、以及在 `q_i = 0` 处有正功率的情形）。
