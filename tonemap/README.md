# tonemap

曝光一致的色调映射，`make test` 跑自检。

## 流程

1. **线性曝光缩放**：`c' = c * exposure`。`exposure` 是线性增益（按档调节时传 `2^ev`），
   只做乘法，通道间比例不变。`exposure <= 0`、NaN、Inf 一律拒绝，函数返回 -1 并清零输出。
2. **只压亮度、保持色相**：取 Rec.709 亮度 `L = 0.2126r + 0.7152g + 0.0722b`，
   对 L 做压缩后得到单一像素增益 `g = curve(L) / L`，三通道同乘 g。
   输出与曝光后颜色严格共线，色相角只残留浮点级误差。
3. **解析曲线**：Reinhard `L/(1+L)`，double 精度版作为参考实现。

## 指标

- `curve_err_abs_max` / `curve_err_rel_max`：float 流水线相对 double 解析参考的
  逐点最大绝对 / 相对误差，自检要求分别小于 `1e-5` 上界。
- `hue_shift_rad_max`：压缩前后归一化 RGB 方向的夹角（弧度），自检要求小于 `1e-5`。
