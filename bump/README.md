# bump

`make test` 跑自检。

高度场 `z = h(x, y)` 的自遮挡（凹凸自阴影）。着色点在原点，`h(0,0) = 0`，patch 为 `[-1,1]^2`、边界上平滑收到 0。沿出射方向 `L`（上半球）步进，比较 `z(t) = t` 与 `h`：

- **按斜率定步长**：设 `lam = max|grad h|`，一步内高度最多变化 `lam * |Lh| * dt`，给定高度预算 `eps` 得 `dt = eps / (lam * |Lh|)`；按步数说就是 `n = ceil(t_exit / dt)`。`|Lh|` 是射线每单位高度对应的横向位移。
- **按步数给精度**：`sweep_steps(id, L, n)` 固定步数；`sweep_slope(id, L, eps)` 按斜率换算步数。步数加倍，误差上界随之缩小。
- **与几何阴影比对**：同一个高度场再用参考光线追踪（固定细步长 + 每个变号区间二分）求真实遮挡，逐方向比较可见性。sweep 只在采样点严格低于表面时报遮挡，因此不会凭空造阴影；误差只可能是漏遮挡。
- **误差上界**：漏掉一个遮挡区间，要求整段遮挡弦都落在两步之间（`maxchord < dt`）。测试按步数 `4 → 1024` 倍增打印 `missed`（实际漏判）和 `bound`（弦短于 `dt`、可能被漏判的阴影射线数）——两者单调下降，`n = 32` 时均为 0；且所有漏判都必然满足 `maxchord < dt`。
- **空高度场拒绝**：`field_ok()` 对不存在的场（负 id、越界 id）返回 0 并写出原因。

## 自检项

| 测试 | 含义 |
|---|---|
| `fields_ok` | 高度场可用，陡峭场的斜率确实更大 |
| `empty_field_rejected` / `empty_field_reason` | 空场拒绝且给出原因 |
| `flat_never_shadows` | 平面在任何方向、任何步数下都不自遮挡 |
| `up_clear` | 正上方永远可见 |
| `reference_stable` | 参考追踪在 2048 / 8192 两步长下结论一致 |
| `scene_has_shadows` | 场景确实有自阴影（否则其余测试无意义） |
| `sweep_never_invents_shadow` | 不凭空造阴影 |
| `miss_implies_short_chord` | 漏判必然对应短于一步的遮挡弦 |
| `error_bound_shrinks` | 步数加倍，漏判数与可解释数单调下降，最终为 0 |
| `eps_halves_steps_increase` | 高度预算减半，步数不减少 |
| `slope_matches_fixed` | 斜率换算的步数与等价固定步数结果一致 |
| `slope_budget_ok` | `eps = 0.01` 下漏判极少 |
