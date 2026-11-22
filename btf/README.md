# btf

按「表面高度 × 视角」双向纹理参数化。`make test` 跑自检（32 项）。

高光方向错在：贴图直接用各向同性的 (u,v) 采样，没有把高度和视角对观感的拉伸
编进参数里，并且纹理缝两侧一阶导突变，越过 seam 处高光条会翻向。这里用一个
可分离、单调、边界一阶导为零的 warp 取代裸查找。

## 参数化

归一化高度 `h ∈ [0,1]`、归一化视角 `w ∈ [0,1]`：

```
u = a(h) = 6h^5 - 15h^4 + 10h^3        (smootherstep)
v = b(w) = (1 - cos(pi w)) / 2         (sin^2(pi w/2))
```

两者都是 C2、在 (0,1) 严格递增、在两端一阶导为零。

## 可逆性

`a`, `b` 在 [0,1] 上连续、严格递增，值域分别为 `[0,1]`，故各自是双射（自带
反函数，代码用二分求逆，往返误差 < 1e-12，`inverse_roundtrip` 逐点验证）。
2D 映射是两条双射的直积，因此 `[0,1]^2 → [0,1]^2` 是双射；`map_u_matches` /
`map_v_matches` 与闭式逐点吻合。

## 雅可比与奇异条件

```
J = [[a'(h), 0], [0, b'(w)]] ,   det J = a'(h) b'(w)

a'(h) = 30 h^2 (1-h)^2            仅在 h = 0, 1 为零
b'(w) = (pi/2) sin(pi w)          仅在 w = 0, 1 为零
```

- 开方内 `h,w ∈ (0,1)`：`det J > 0`，2×2 雅可比非奇异，映射局部可逆
  （`jacobian_interior_positive`）。
- 奇异集恰为 `[0,1]^2` 的边界（`jacobian_boundary_zero`）；四角两轴同时奇异。
- 2-范数条件数 `kappa = max(|a'|,|b'|) / min(|a'|,|b'|)`。中心处
  `a'(1/2)=1.875`、`b'(1/2)=pi/2`，`kappa ≈ 1.193`（`cond_center_value`）。

## 边界 C1

`a'(0)=a'(1)=0`、`b'(0)=b'(1)=0`（`a_boundary_c1` / `b_boundary_c1`，按容差
判零，因为 `sin(pi)` 有 1e-16 量级浮点残差）。于是导数在 seam 两侧都取 0，
一阶连续 —— 跨 seam 高光条方向不再跳变。二阶导在边界有限
（`boundary_second_deriv_finite`），故整体 C2、边界 C1 成立。

## 与逐点采样比对

把每个格子的四角按精确映射取出、做线性插值（点阵查找在格内就是这样重建），
测它与精确映射的偏差。因映射可分离、每个输出只依赖一个轴，2D 双线性误差退化为
该轴的一维线性插值误差。上界用拉格朗日余项

```
|f - L| <= (step^2 / 8) * max|f''|
max|f''| <= |f''(格首)| + step * max|f'''|
max|a'''| = 60,   max|b'''| = pi^3 / 2
```

逐个格取 max 得到严格上界 `sample_error_bound_holds`（实测 ≤ 上界），并验证
细化加密时误差二阶收敛（`sample_error_second_order`，32→64 误差比 > 3.5）。

## 条件数上限与近奇异报告

`KAPPA_MAX = 10`。`cond_report(n)` 扫 n×n 网格：

- 边界点 → 记入 `n_boundary_singular`（精确奇异）；
- 内点 `kappa > KAPPA_MAX` → 记入 `n_near_singular`，并记录最坏处
  `(h_worst, w_worst)` 与 `kappa_worst_near`；
- 其余 → `n_admissible`，`kappa_max_admissible` 必须 `<= KAPPA_MAX`。

自检确认上限被守住、且近奇异点确实被报告（`report_ceiling_respected` /
`report_near_singular_found` / `report_worst_near_above_ceiling`）。

## 零维拒绝

`param_init` 拒绝 `dims <= 0`（零维）、`dims != 2`、以及任一区间被钉成一点
（高度或视角范围坍缩，映射退化为 1D/0D）。被拒的句柄 `param_map` 返回 0，
不会悄悄产出 1D 查找。对应 `reject_dims_zero` / `reject_dims_one` /
`reject_constant_height` / `reject_constant_view` / `rejected_handle_no_map`。
