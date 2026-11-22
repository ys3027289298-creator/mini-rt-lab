# material

`make test` 跑自检。

## 互易性校验

`make test` 会为每个材质跑一次互易性扫描（4096 对方向）：

- 方向对由 BRDF 自身采样产生：先绕法线采 `wi`，再以 `wi` 为条件采样 `wo`。
- 对每对方向分别求值 `f(wi, wo)` 与 `f(wo, wi)`，统计绝对偏差与相对偏差
  （`2|a-b|/(|a|+|b|)`）的分布：exact、≤1e-6、≤1e-3、≤1e-1、>1e-1，
  以及 p50/p90/p99 与最大偏差。
- `class=reciprocal_analytic` 的材质**必须**得到 `max_abs_dev=0`、全部落在
  exact 桶（`verdict=exact_zero`）。这不是容差判定：求值按对换对称形式书写
  （共享 `wi+wo`、对称的 `G`/`1-|wi·h|`/分母），两个参数顺序走完全相同的
  浮点运算，所以要求的是逐位相同。任何非零偏差判 `FAIL_nonzero_deviation`。
- `class=nonreciprocal` 的材质单独标注为 `labeled_nonreciprocal`，不算错误。
  只要扫描在“有信号”的方向对上测到偏差即算识别成功。
  `labeled_nonreciprocal_undetected` 表示扫描没测出预期偏差，才需要看。
- `selftest_directional_defect` 是故意做坏的探针：它声明为解析互易却带有
  方向项，输出**应当**是 `FAIL_nonzero_deviation`，自检断言其返回失败，
  用来证明本校验不是空转。

方向向量未归一（容差 1e-4，零向量也算）时不求值，直接报出出错的方向：

```
direction_not_unit material=ggx_r0.35_analytic status=wi_not_unit direction=wi vector=(...) length=2.000000 tolerance=1e-04
```

`status` 取值 `wi_not_unit` / `wo_not_unit` / `both_not_unit`，两个方向都不合法
时分两行分别报出各自的向量和长度。
