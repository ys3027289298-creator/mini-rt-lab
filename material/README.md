# material

`make test` 跑自检。

## 多层材质合成

`Layer { albedo[3], coverage }`，`layers_composite(Layer*, count, LayerReport*)`：

- 输出按各层 **coverage 加权平均 albedo**（不是逐层累加，避免叠层后变亮变白）。
- 各通道结果 clamp 到 `[0,1]`；总反射率不超过 1（能量不增）。
- **覆盖率之和 > 1** 时按比例归一并记录：`coverage_sum` 存原始和，`norm_scale`
  存 `1/coverage_sum`，`normalized=1`；`layer_cov[]` 存归一后实际权重。
- 输入校验：coverage 为负 / 非有限、albedo 越界 `[0,1]`、NaN、层数为 0 或超过
  `LAYER_MAX`(=8) 均返回 0。

对比与误差上界：

- `layers_reference`：逐层解析 "over" 合成，第 i 层权重
  `cov_i · Π_{j<i}(1-cov_j)`，不做归一。
- `layers_error_bound`：给出合成值与该参考值之差的逐通道上界
  （由极值 albedo 与正/负 coverage 偏差累积得到，含 float 舍入余量）。

`artifacts/layer_probe.c` 可直接编进本文件做端到端验证：

```sh
gcc -O2 -std=c11 layer_probe.c -o probe -lm && ./probe blend|partial|energy|bound|invalid
```
