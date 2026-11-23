# mini-rt-lab

软件光线追踪基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `bvh/` | 包围盒树构建、中位数分裂、最近命中 |
| `sample/` | 分层与低差异采样、方差估计 |
| `material/` | Lambert 与 GGX 求值、余弦加权采样 |
| `shared/` | 统一随机源（PCG32）与三模块共享场景 |

## 统一随机源

三个模块原先各用各的随机/准随机机制（Halton、网格算术、中位数排序），
换机器或换编译器后行为无法统一对照。现在所有随机来源收敛到
`shared/rng.h` / `rng.c` 的 PCG-XSH-RR 64/32：

- 定宽整数运算（`stdint.h`），不依赖 libc `rand()`、时钟或平台状态，
  同一 `(seed, stream)` 在任何机器/编译器上产生逐位相同的序列。
- `shared/scene.h` 定义唯一共享场景：一个 `SCENE_ID`、一个 `SCENE_SEED`，
  采样/材质/几何分别使用 stream 1/2/3，序列互不重叠且同源。

编译开关 `RT_UNIFIED_RNG`（模块默认开启）。关闭后每个模块走原始代码路径，
不引入任何额外输出，行为与改动前**逐字节相同**（见下 `legacy-check`）。
公开函数签名均未改动；统一路径只在原有函数内部接入 RNG。

## 三处验收（缺一份数据即未完成）

同一场景（同 seed、不同 stream）驱动三份数据，全部落在
`artifacts/acceptance/`：

| 模块 | 验收方式 | 产物 |
|---|---|---|
| `sample/` | 均匀性直方图 + 卡方/相对偏差 | `histogram.csv` |
| `material/` | cosθ 样本方差对解析值 1/18 | `variance.json` |
| `bvh/` | 场景渲染图 + FNV-1a64 字节哈希 | `image.ppm`、`image_hash.json` |

每个文件头部都带 `scene_id`、`seed`、`stream`，由 `verify_artifacts.py`
核对三份数据确属同一场景、三条独立 stream，并独立复核直方图计数、
方差容差与图像哈希（用 Python 重算 FNV-1a）。

## 常用命令

```sh
make verify         # 生成三份数据并交叉核对（RESULT: PASS）
make legacy-check   # 关闭统一来源，三模块 stdout 与 golden 逐字节 cmp
make cross-check    # clang 重新构建，四份产物与 gcc 逐字节相同
make clean

# 单模块
make -C sample test          # unified 自检
make -C sample test-legacy   # 原始代码路径自检
```

失败语义：`verify_artifacts.py` 在缺文件、场景/seed 不一致、stream 不是
1/2/3、直方图不均匀、方差超差、哈希对不上时均以非零码退出。
