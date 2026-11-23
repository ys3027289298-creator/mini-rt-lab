# denoise

采样器 ↔ 降噪器统一帧接口、时间累积与闪烁自检。`make test` 跑自包含自检。

## 背景：三处接口错位导致画面忽明忽暗

旧实现双方各自声明结构、各自解释字段，没有共享契约：

| # | 错位点 | 采样器侧 | 降噪器侧 | 后果 |
|---|---|---|---|---|
| 1 | 载荷语义 | 字段名叫 `sum`，实际填的是**均值** | 当样本和再除一次 spp | spp=4 时亮度被压到 1/4、spp=8 时 1/8，帧间明暗交替 |
| 2 | 方差语义 | 发的是**单样本方差** | 当成均值方差直接用于混合权重 | 权重不随样本数收缩，噪声永不收敛 |
| 3 | 帧计数 | 0-based | 缓存计数从 1 开始，按自己的计数反推 spp | 奇偶 spp 被取反，错位叠加 |

spp 在 4/8 间交替时，旧路径输出随帧周期性明暗振荡；自检复现到 **27/27 个有效帧间间隔闪烁**。

## 统一接口（`src/rt_frame.h`，唯一契约来源）

### 接口兼容清单

| 符号 | 方向 | 类型/单位 | 统一语义 | 旧路径对应字段 | 是否兼容旧语义 |
|---|---|---|---|---|---|
| `RtFrame.mean` | 采样器→降噪器 | `const float*`，线性辐射度均值 | 已按本帧 spp 归一化的**均值** | `LegacySamplerPacket.sum`（名实不符） | 否，禁止再除 spp |
| `RtFrame.var` | 采样器→降噪器 | `const float*`，单样本方差 | Bessel 修正，与 `sample/` 模块 `Stat.var` 同公式 | `varmean`（被误当均值方差） | 否，均值方差由降噪器推导 |
| `RtFrame.spp` | 采样器→降噪器 | `int`，本帧样本数，>0 | 降噪器只信该字段，不自行反推 | 无（降噪器用私有帧计数猜） | 否，消除猜测 |
| `RtFrame.frame` | 采样器→降噪器 | `int`，0-based | 唯一帧序号 | 双方各自计数 | 否，单一计数源 |
| `RtFrame.width/height` | 采样器→降噪器 | `int` 像素 | 行优先、stride==width、无填充 | 双方各自缓存 | 形状不匹配 `rt_accum_push` 返回 -1 |
| `RtAccum.mean` | 降噪器内部/输出 | `float*` | 按总样本数 N 加权的累计均值 | 无 | 新增 |
| `RtAccum.vpool` | 降噪器内部 | `float*` | `sum(spp*var)` 单样本方差池 | 无 | 新增 |
| `RtAccum.mvar` | 降噪器内部/输出 | `float*` | 累计**均值**方差 = `vpool / N²` | 无 | 新增（关键修正：不是 `/N`） |
| `RtAccum.samples` | 降噪器内部 | `long`，累计样本总数 N | 跨帧唯一样本计数 | 无 | 新增 |
| `rt_accum_push` | 调用 | 返回 0/-1 | 形状校验 + 样本数加权递推 | `legacy_denoise` 固定权重混合 | 新契约，旧函数仅留作复现 |
| `rt_frame_passthrough` | 调用 | `memcpy` | 关闭统一时逐像素原样输出 | 无 | 旁路，字节级一致 |

兼容规则：

1. 双方只允许 `#include "rt_frame.h"`，禁止再声明同名私有结构。
2. 采样器永远发均值和单样本方差；任何 `/spp`、`/N` 只能出现在降噪器内。
3. `rt_accum_push` 形状不匹配返回 -1，不静默截断。
4. 关闭统一走 `rt_frame_passthrough`，输出与输入逐字节一致（有 f32 产物与独立脚本核验）。

## 复现步骤

```bash
cd denoise
make clean && make test          # 自包含：test_out/ 内生成 golden 再校验，all_ok=1

# 交付产物（固定目录，确定性输出，跨构建字节一致）
make
DENOISE_OUT=../artifacts/denoise ./denoise --write-golden   # 首次：冻结旧路径 golden
DENOISE_OUT=../artifacts/denoise ./denoise                  # 之后：对照 golden 校验

# 独立核验：结论从产物数据重算，不信任 C 自检输出
python3 -I tools/verify_artifacts.py ../artifacts/denoise   # OVERALL=PASS

# PGM -> PNG 便于肉眼检查
python3 -I tools/pgm_to_png.py ../artifacts/denoise/noisy_frame00.pgm /tmp/noisy.png
```

## 产物清单（`artifacts/denoise/`）

| 文件 | 内容 |
|---|---|
| `noisy_frame00.pgm/.png` | 噪点图像摘要：第 0 帧原始 MC 输出（64×48，双高斯光斑 + 高斯噪声） |
| `legacy_output_frame39.pgm/.png` | 旧路径末帧：整体偏暗、残留噪点 |
| `unified_output_frame39.pgm/.png` | 统一路径末帧：平滑、亮度正确 |
| `variance_stats.csv` | 40 帧逐帧方差统计：两路整图均值、对真值 MSE、单样本方差、统一路径均值方差 |
| `golden_legacy_off.f32` | 旧路径末帧冻结输出（magic `RTG1` + w/h + float 载荷），关闭统一的字节级基线 |
| `bypass_input.f32` / `bypass_output.f32` | 旁路输入/输出，载荷逐字节一致 |

## 关键结论（自检实测）

- 旧路径闪烁帧 **27**，统一路径闪烁帧 **0**（暖机 12 帧，阈值 2e-3 相对变化）。
- 旧路径偶/奇帧亮度比 **1.1001**（EMA 阻尼后的残余振荡），统一路径 **1.0000**。
- 统一路径均值方差从 0.002440 降到 0.000081（约 1/30，符合 N 增长），对真值 MSE 0.000110 vs 旧路径 0.0407。
- 旁路输入/输出 3072 个 float 全部字节一致；golden 跨独立构建 md5 一致。
- 确定性来源：xorshift64 固定种子，两路共享同一样本序列，差异只来自接口处理。
