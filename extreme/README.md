# extreme/ — 极端场景矩阵

对 `bvh/`、`sample/`、`material/` 三个模块做极端场景压力测试。
**不改动仓库样例**：驱动通过 `#include` 直接引用模块源码（`main` 重命名），
构建产物只落在 `extreme/build/`；`make verify-clean` 校验已跟踪文件零改动。

## 方法

- 每个用例在独立子进程运行（`harness.h`）：崩溃（信号 / sanitizer 中止 /
  10s 挂起）被隔离计数，不中断矩阵。
- 编译带 `-fsanitize=address,undefined -fno-sanitize-recover=all`，
  UB 当场转为可计数的崩溃。
- 所有输出经 `EXPECT_FINITE` 检查 NaN / Inf，经 `EXPECT_TRUE` 比对人工构造的
  期望值（baseline 用例含手算精确值，如 `tmin=1.0`、GGX NDF=5.0930）。
- 判定标准：**crashes == 0 且 violations == 0** 才算通过。

## 运行

```sh
make check         # 构建并跑全部 34 个用例，任一崩溃/NaN/Inf 则退出码非零
make verify-clean  # 确认仓库已跟踪文件零改动
```

## 极端场景类别与触发条件

| 类别 | 触发条件 | 用例 |
|---|---|---|
| 超长细三角形 | 包围盒边长比 ≥ 1e12；或薄边尺寸低于坐标量级的 float 相对精度（~1e-7，如 1e-6 的边在 1e6 坐标下被舍入吞掉） | `long_thin/*` |
| 极近/远平面 | tmax ≤ 1e-30（极近）或 ≥ 1e30（极远）；tmax=0；射线原点距图元 ≥ 1e25 | `near_far/*` |
| 零面积光源 | 包围盒 min==max（零体积/零面积）；GGX roughness=0（delta 分布 = 点光源极限）；rough⁴ < FLT_EPSILON/2（rough < ~0.016，d 坍塌除零）；掠射角 cos=0（pdf=0）；采样数 strata=0 | `zero_area/*` |
| 超大坐标 | 坐标 ≥ 1e18（float 精度坍塌）；坐标 ≥ 1.7e38（质心 min+max 溢出）；包围盒边长 ≥ 1e19（面积溢出）；被积函数 ≥ 1e308（double 累加溢出）；序列索引取 UINT_MAX | `huge_coords/*` |

## 实测结果（de9c470，gcc -O1 + ASan/UBSan）

**崩溃：0 / 34。** 所有类别（含 32 个同质心点图元建树的无限递归风险、
次法向方向分量、tmax=0）均正常返回，无信号、无 sanitizer 报告。

**NaN/Inf：9 处违规，分布在 7 个用例** —— 按"崩溃与 NaN 必须为零"的标准，
当前代码**不通过**，需修复或文档化以下缺陷：

| 模块 | 用例 | 现象 | 触发条件 |
|---|---|---|---|
| material | `ggx_rough_zero` | NaN（0/0） | roughness = 0 |
| material | `ggx_rough_0p01` | Inf | rough⁴ < FLT_EPSILON（rough ≲ 0.016，**含常用粗糙度 0.01**）：`1-(1-a2)` 舍入为 0，d 坍塌除零 |
| material | `ggx_rough_underflow_nan` | NaN | rough < ~1e-19：rough⁴ 本身下溢为 0 |
| bvh | `near_fltmax` | Inf ×3 | 坐标 ≥ ~1.7e38：质心 `(min+max)*0.5` 溢出 |
| bvh | `area_overflow` | Inf | 边长 ≥ ~1e19：`dx*dy` 溢出 |
| sample | `zero_strata` | NaN | strata = 0：`sum/n` 即 0/0 |
| sample | `exploding_integrand` | Inf | 被积函数 ≥ 1e308：累加溢出（var 的 NaN 被 `var>0` 守卫吞掉，mean 仍漏出 Inf） |

## 已验证安全的边界（人工构造比对通过）

- 边长比 1e12 的超长细盒：建树、命中、薄边外 miss 全部正确；
  薄边被 float 舍入吞掉后退化为零厚度，仍不崩。
- tmax=0 / 1e-30 / 1e30、原点在盒面上、原点距 1e25：hit/miss 与 tmin 值精确。
- 零体积点盒、零厚度平板：命中与 tmin 精确，面积为 0 不 NaN。
- 32 个质心完全相同的图元：分裂落空走 count/2 兜底，建树终止且节点数 ≤ 64。
- 原点 -3e38 对盒 +3e38：跨度溢出为 inf 后 slab 测试凑巧给出正确 miss。
- `cosine_sample` 边界 (u∈{0,1}) 输出单位向量；u1<0 越界时 NaN 被
  `norm3` 的 `l>0` 守卫吞掉返回零向量（安全但静默，值得注意）。
- `cosine_pdf(0)=0`：本身有限，但下游 1/pdf 即 Inf —— 调用方需自护。
- Halton 索引取 UINT_MAX（含 radical_inverse 溢出回绕路径）：输出有限且在 [0,1)。
