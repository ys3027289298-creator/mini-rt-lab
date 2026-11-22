# invariants — 渲染器物理不变量检查

对仓库里的三个渲染模块（`bvh/`、`sample/`、`material/`）做**只读**的物理不变量审计。
检查器不修改模块源码、不调用模块自己的 `check()`/`main()`、不写模块状态；它把三个
源文件**原样**编进观察用的 shim（仅重命名入口点），然后只调用其中的纯函数。

```
make validate   # 先证明“尺子是准的”：无副作用 + 无悔报 + 违例必检出
make audit      # 再用尺子量渲染器（退出码 1 = 发现违例，这是结论而非故障）
make check      # 两者按顺序跑
```

`make` 下任何非零都显示为退出码 2（make 的约定）；直接运行二进制能得到真实退出码。
直接运行：`./invcheck`、`./invcheck --calibrate`、`./invcheck --render-check`。

## 四类不变量与判定阈值

所有阈值都是**判定阈值**，不是拟合值：每个都取在该检查数值噪声底之上至少两个数量级、
在标定用的注入违例之下至少一个数量级，同一阈值要同时把“正确”和“破坏”分开。
阈值集中定义在 `include/invariants.h`，每次运行会把实测值与阈值并排打印，便于复核。

| 类 | 检查 | 阈值 | 含义 |
|---|---|---|---|
| **energy** | `lambert_reflectance_le_one` | 1e-3 | Lambert 反射率 = albedo，要求 ≤ 1 |
| | `lambert_mc_consistency` | 1e-2 | 用模型**自己的**采样器+pdf 估的反射率 = 解析值 |
| | `brdf_nonnegative` | 1e-6 | BRDF 非负 |
| | `cosine_pdf_normalized` | 2e-3 | ∫ 余弦 pdf 立体角 = 1 |
| | `ggx_ndf_normalized` | 5e-3 | 2π∫ D(μ)μdμ = 1（n=2e5 Simpson） |
| | `ggx_finite_nonnegative` | 0 bad | 合法域内 D 不得为 NaN/Inf/负 |
| **geometry** | `box_hit_vs_exact` | 0 | 与精确 double 判交比较，清晰区不得有分歧 |
| | `boundary_hit_recall` | 假阴性 0 | 清晰真命中不得漏检 |
| | `boundary_false_hit` | 假阳性 0 | 清晰真未命中不得报命中 |
| | `hit_t_accuracy` | 2e-3 | 命中距离相对误差 |
| | `tree_*` | 1e-4 | 根箱覆盖、子箱包含、图元恰好引用一次 |
| | `ambiguity_band_exercised` | >0 | 切触用例确实落入不确定带 |
| **sampling** | `sample_domain` | 0 | 采样点必须在 [0,1]，>1 或 <0 为硬违规 |
| | `uniform_mean` / `uniform_variance` | 1e-2 / 5e-3 | 均值 1/2、方差 1/12 |
| | `discrepancy` | 1e-2 | 一维星偏差 D* |
| | `distinctness` / `integral_consistency` | 1e-3 / 2e-2 | 去重率、MC 积分 ∫u²=1/3 |
| **symmetry** | `sampler_azimuthal_mean` | 1e-2 | E[x]=E[y]=0（方位对称） |
| | `sampler_var_equal` | 2e-2 | Var[x]=Var[y] |
| | `sampler_mirror_pointwise` | 1e-5 | s(u1,u2) = reflect_xz s(u1,1−u2) |
| | `ggx_even_in_nh` | 1e-6 | D(h)=D(−h)（只依赖 \|n·h\|） |

关于反射对称性：本仓库唯一的 BRDF 是常数 Lambert，**Helmholtz 互易性在这里是空条件**，
所以对称性检查落在唯一有内容的地方——采样器的方位对称/镜像恒等式，以及 GGX 关于 n·h 的偶性。

几何类的“不确定带”：float 判交在面/边/切触附近本就可能与精确算术不同。检查器对每个
用例做单 ulp 扰动，判定其答案是否稳定；**只在稳定（清晰）区**要求与精确算术一致。
这样阈值不会把合法的 float 舍入误判为缺陷，同时仍能抓出真正的逻辑错误。

## 标定（“违例必被检出且误报为零”）

`--calibrate` 用**人工构造的违例**（`src/mutants.c`，每个只在黄金参考上破坏一处）验证尺子本身：

- 10 个违例全部被**目标类**检出（`false_negatives=0/10`）；
- 两个已知正确实现（double 黄金参考 + **纯 float 控制实现**）**零误报**；
- 每个违例不牵连其它类（`collateral=0`）。

纯 float 控制（`inv_float_control`）是关键：它与被审计模块同为 float 精度、但实现正确。
若它失败，则任何归咎于渲染器的失败都不可信。它能通过，说明阈值容得下正常 float 误差。

> 标定过程中发现过一个**检查器以外**的问题并已修正：最初用 `1−halton(i)` 作为采样违例，
> 但它在 [0,1) 上是保测度的对称变换（均值/方差/偏差/积分不变），**不是** 1 维违例——
> 真正的 1 维违例是重复（`halton(i/2)`）或被压缩（`0.5·halton`）。已换成真正的违例。

## 无副作用证明

`--render-check` 对一段固定“渲染工作负载”（建树 + 射线扇 + 采样 + BRDF 求值）做 64 位
FNV 哈希，在跑完整个不变量检查**前后**各取一次，要求逐位相同；并读取三个模块各自的
`failures` 计数器，要求在一次检查后仍为 0。结果：

```
bit_identical=1
module_state_untouched=1
verdict=NO_SIDE_EFFECTS
```

## 审计结论（当前 HEAD `de9c470`）

尺子先通过（黄金参考、float 控制均 clean，无副作用），然后测得渲染器存在真实违例：

| 类 | halton2 | halton | 结论 |
|---|---|---|---|
| energy | **违例** | **违例** | GGX 在 `rough=0` 且 `nh=1` 返回 NaN；小 `rough` 返回 Inf |
| geometry | **违例** | **违例** | 近切射线被漏检（见下） |
| sampling | ok | **违例** | `halton()` 高位索引越界 |
| symmetry | ok | ok | — |

具体根因（均已用独立探针复现）：

1. **几何漏检**（9/20149 清晰用例，全部在 `boundary_hit_recall`）：
   `bvh/src/bvh.c:39` 的 `fabsf(d[i]) < 1e-8` 把很小的方向分量当作“平行”，随后只在
   原点位于 slab 外时判 miss。对**从很远处以极小角度命中盒子**的射线（如
   `dist=1e9`, `d_y≈−1.2e-9`），原点在 y-slab 外、但射线确实穿过盒子——于是被误判漏检。
   这是边界采样命中率直接失守。合法 float 实现（float 控制）零漏检，证明可修。

2. **能量：GGX 非有限**：`ggx_ndf(1.0, 0.0)` 得 `-nan`（0/0），`rough∈{1e-4,1e-3,0.01}`
   在 `nh=1` 得 `inf`。渲染器把这些值累加进辐射度会污染整条路径。正确实现应在奇异点取有限约定（如 0）。

3. **采样越界**：`sample/src/sample.c:11-17` 的 `radical_inverse` 用 `unsigned` 累乘
   会溢出，`halton()` 又额外 `× base / 2^32`。高位索引出现 `1.39、2.09` 等**超出 [0,1)**
   的值（7 个）。注意模块自带测试只测了 `halton2()`，从未测过 `halton()`，后者漏网。

Lambert/余弦 pdf、树结构包含与划分、采样器对称性与 GGX 偶性均通过。

## 布局

```
invariants/
├── Makefile                 # validate / audit / check
├── include/invariants.h     # 公共 API + 全部阈值
├── include/inv_impl.h       # 内部：黄金参考 / shim 适配器
└── src/
    ├── shim/edge_bvh.c      # 原样 #include bvh/src/bvh.c（仅改入口名）
    ├── shim/edge_sample.c   # 原样 #include sample/src/sample.c
    ├── shim/edge_material.c # 原样 #include material/src/material.c
    ├── golden.c             # double 黄金参考 + 纯 float 控制
    ├── mutants.c            # 10 个人工违例（标定用）
    ├── inv_checks.c         # 四类不变量的实现
    ├── renderer_edge.c      # 把三个 shim 组装成 InvModel
    └── inv_main.c           # 驱动：audit / calibrate / render-check
```
