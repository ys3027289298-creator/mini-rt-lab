# alloc

按像素复杂度分配维度（每像素样本数），`make test` 跑自检。

## 问题

所有像素拿同样多维度时，低复杂度像素的估计早就收敛、维度被浪费。
本模块在**总维度预算 B 固定**的前提下，把维度向高复杂度像素倾斜。

## 模型

- **复杂度**：试点采样（pilot）的无偏样本标准差 `pilot_std()`，
  `s = sqrt(sum (x - mean)^2 / (k-1))`，k<2 或常量场景返回 0。
- **误差与维度关系**：像素 i 拿 n_i 个维度时，均值估计的 MSE
  `v_i(n_i) = sigma_i^2 / n_i`（独立样本、方差不变时）。
- **优化目标**：最小化总 MSE `E = sum sigma_i^2 / n_i`，
  约束 `sum n_i <= B` 且每个 `n_i >= n_min`。
- **解**：拉格朗日求导得 `n_i ∝ sigma_i`（水位法 water-filling，
  `n_i = max(n_min, sigma_i / lambda)`），取整后按边际收益
  `sigma_i^2 / (n_i(n_i+1))` 贪心发剩余维度。

## 与均匀维度对比

- 均匀参考：`E_uni = (N/B) * sum sigma_i^2`
- 连续最优比例（柯西不等式上界）：
  `E_opt / E_uni = (sum sigma_i)^2 / (N * sum sigma_i^2) <= 1`
  复杂度越不均匀收益越大；全均匀场景比例恰为 1。
- 自检会打印逐像素分配表、两种 MSE、比例与下降百分比，并扫多个预算验证。

## 拒绝规则（维度为零一律拒绝）

- `budget < npix * n_min`（含预算为 0）：无法保证每像素至少 n_min 维 → `ALLOC_EBAD`
- `n_min < 1`、`npix` 非法、`sigma` 为负或 NaN → `ALLOC_EBAD`
- 成功时断言：每个 `out[i] >= n_min`，且 `sum(out) == budget`
  （不超预算、预算用尽、无零维像素）

零复杂度像素（sigma=0）合法：只拿 n_min 个维度，其预测 MSE 为 0。
