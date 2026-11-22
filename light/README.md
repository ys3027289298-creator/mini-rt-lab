# light

锥体光源采样与距离立体角修正，`make test` 跑自检。

## 模型

半径 `R` 的圆盘面光源，从着色点（距离 `d`）看张开**精确**立体角

```
Omega_disk(d) = 2*pi*(1 - 1/sqrt(1 + (R/d)^2))
```

它恰好等于半角 `alpha_exact = atan(R/d)` 的锥体立体角：

```
Omega_disk(d) == Omega_cone(atan(R/d))
```

因此锥体近似在参考距离处无误差。锥体在 `d_ref` 处按 `alpha_ref = atan(R/d_ref)`
拟合，距离修正按 `1/d^2`：

```
Omega_corr(d) = Omega_cone(alpha_ref) * (d_ref/d)^2
```

## 采样

`cone_sample` 在修正后的有效锥体内**均匀按立体角**采样：

- `cos(theta)` 在 `[cos(alpha_eff), 1]` 上均匀，`phi` 在 `[0, 2*pi)` 上均匀；
- pdf（每立体角）为 `1 / Omega_cone(alpha_eff)`；
- 半角为零时退化为点（delta）光源：方向取轴，pdf 为 `+inf`。

`alpha_eff` 由修正后的立体角反解：`cos(alpha_eff) = 1 - Omega_corr/(2*pi)`。
近距时 `1/d^2` 会过冲；一旦修正立体角达到 `2*pi`，`cone_light_alpha_eff`
返回 `ERR_SOLID_ANGLE`。

## 与精确面积采样比对

`cone_width_compare` 给出圆盘精确半角 `alpha_exact`、修正锥体有效半角
`alpha_eff`，以及软影宽度比

```
ratio = tan(alpha_eff) / tan(alpha_exact),   err = |ratio - 1|
```

宽度取半角正切，是因为直边遮挡下的半影半宽与 `tan(alpha_light)` 成正比
（系数只取决于遮挡几何，与光源模型无关，故比值即模型误差）。

`disk_solid_angle_area_mc` 用面积均匀采样对圆盘立体角做 Monte Carlo
交叉校验（自检里与闭式解比对，误差 < 1%）。

## 软影宽度误差上界

设 `t = R/d`、`q = Omega_corr/Omega_disk = F(t_ref^2)/F(t^2)`，
其中 `F(s) = 1/(sqrt(1+s)(sqrt(1+s)+1))`。宽度比满足

```
tan^2(alpha_eff)/tan^2(alpha_exact) - 1 = (q - 1) * A(q, u),
A(q, u) = (2 - (q+1)u) / ((1 - q u)^2 (2 - u)),   u = 1 - cos(alpha_exact).
```

在 `alpha_cap` 内取 `|q-1| <= Bq`、`|A| <= Amax`，得**先验上界**

```
err <= Bw = 1 - sqrt(1 - Bq*Amax).
```

`cone_width_error_bound(alpha_cap)` 返回 `Bw`。有效范围 `alpha_cap <= 35°`
（闭式解在约 40.3° 附近进入不合适区间；这里留出裕度）。自检在所有
参考角与距离上验证 `err <= Bw`、`Bw` 随 `alpha_cap` 单调递增。示例：
`alpha_cap = 30°` 时 `Bw ≈ 0.21`。

## 越界报告

`LightStatus`：

| 状态 | 含义 |
|---|---|
| `ERR_ANGLE_RANGE` | 半角 `< 0`、`>= pi/2`，或 `alpha_ref != atan(R/d_ref)` |
| `ERR_ANGLE_BOUND` | 半角在物理范围内但超出上界有效域（`> 35°`） |
| `ERR_SOLID_ANGLE` | 距离修正后立体角 `<= 0` 或 `>= 2*pi` |
| `ERR_UV` | 采样参数越界 |
| `ERR_AXIS` / `ERR_RADIUS` / `ERR_DISTANCE` / `ERR_NULL` | 参数非法 |
