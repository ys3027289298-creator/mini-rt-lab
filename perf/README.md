# perf -- 三处渲染路径的统一性能预算、两档验收

三处渲染路径（`bvh/` 几何、`sample/` 采样、`material/` 材质）各自有预算，
但从未被量过。这个模块把预算收敛到一份文件 `budget.conf`，用**低/高两档**
对同一套路径做验收：

| 验收项 | 低档 | 高档 | 判定 |
|---|---|---|---|
| 分辨率 | 80×60 | 160×120 | — |
| 每像素 rays(BVH) / draws(sample) / samples(material) | 4 / 64 / 8 | 16 / 256 / 32 | — |
| 工作单位（= w×h×spp） | 19.2k / 307.2k / 38.4k | 307.2k / 4.92M / 614.4k | 必须 **等于** 上界 |
| 单路径耗时中位数 | ≤ 100 ms | ≤ 400 ms | 不得超上界 |
| 全档耗时中位数 | ≤ 250 ms | ≤ 800 ms | 不得超上界 |
| 峰值内存 (RSS) | ≤ 32 MiB | ≤ 64 MiB | 不得超上界 |
| 图像哈希 | 每路径 + 三图合并 | 同左 | 必须可复算 |

**缺任一档数据即判未完成**：`verify_tiers.py` 找不到任一档的
`metrics.json`、任一 `.ppm`，或缺字段，都直接 FAIL，不会把"没跑"当成"通过"。

## 三处路径怎么测

每个 shim 都把对应模块的源文件**原样 include**（只把模块自带的 `main`
改名），所以测的是仓库里的真实代码：

- `perf_bvh.c`：`bvh/src/bvh.c` 的 `tree_add` / `build` / `aabb_hit`，
  每像素 `spp` 条抖动光线做最近命中；
- `perf_sample.c`：`sample/src/sample.c` 的 `stratify()` / `f_quad()`，
  每像素 `spp` 次分层抽样，再用模块自己的 `halton2()` 产生空间纹理；
- `perf_material.c`：`material/src/material.c` 的 `cosine_sample()` /
  `brdf_lambert()` / `ggx_ndf()`，每像素 `spp` 个材质样本。

三个 shim 各有一个工作量计数器，每次渲染前清零、每处理一个工作单位
（一条光线/一次抽样/一个样本）加一。`verify_tiers.py` 要求
`measured == scheduled == w×h×spp`，所以"少算了工作"也会被抓出来。

随机性来自 `include/perf_rng.h` 的 PCG-XSH-RR 64/32，只用定宽整数，
同一 `(seed, stream)` 逐位可复现；三条路径用互不相同的 stream 1/2/3。
浮点运算显式 `-ffp-contract=off`，实测 gcc -O0/-O2/-O3 与 clang -O2
生成的 PPM **字节完全相同**。

## 预算检查开关

`perf run <tier> ... --no-budget-check` 只跳过"上界比较"这一步；
渲染参数、种子、像素缓冲与哈希全部不变。`make no-check` 会把关闭检查后
的 PPM 与正常档的 PPM 逐字节比较，要求**逐像素不变**。

## 图像哈希

哈希是该 PPM 文件**全部字节**的 FNV-1a64（含 `P6` 头）。C 端在写盘后
直接读文件计算，`verify_tiers.py` 用 Python 独立重算，任何像素或头部
差异都会让两者不一致。`combined_image_fnv1a64` 是把 bvh/sample/material
三个 PPM 按顺序拼接后的哈希。

## 用法

```sh
make -C perf test        # 自检：RNG 纯度、渲染确定性、工作计数精确
make -C perf low         # 低档 -> artifacts/perf/low/
make -C perf high        # 高档 -> artifacts/perf/high/
make -C perf accept      # selftest + 两档 + 校验 + 关闭检查逐像素不变
make -C perf no-check    # 只做关闭预算检查后的逐像素对比
make -C perf clean
```

产物（`artifacts/perf/`）：

| 文件 | 内容 |
|---|---|
| `low/metrics.json`、`high/metrics.json` | 参数、上界、每路径耗时表、工作单位、峰值 RSS、哈希、`all_ok` |
| `{tier}/{bvh,sample,material}.ppm` | 三条路径的图像 |
| `nocheck-{low,high}/` | 关闭预算检查后的产物（用于逐像素不变证明） |
| `REPORT.md` | 两档耗时表、峰值内存、图像哈希、复现命令（由 verify 生成） |
| `env.txt` | 测量主机、编译器、核数、UTC 时间 |

## 上界怎么来的

工作量上界是确定值（`w×h×spp`），不依赖机器。时间/内存上界取参考机器
（见 `artifacts/perf/env.txt`）实测值的 20–100 倍余量，是**回归闸门**：
它能抓住某条路径变得数量级更慢或分配数量级更多的内存，而不把正常的调度
抖动当成失败。测量取 7 次重复的中位数。

## 失败语义

`perf run`：工作单位与计划不符、任一耗时/内存超界 → 退出码 1；参数或
I/O 错误 → 2。`verify_tiers.py`：任一档缺失、哈希不符、工作/耗时/内存
超界、或关闭检查后像素改变 → 退出码 1。`make accept` 任一步非零即中止。
