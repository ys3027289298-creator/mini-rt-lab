# pipeline

五阶段渲染管线与直通排查工具，依赖仓库内 bvh/material/sample 的同款数学约定。

## 阶段

| # | 阶段 | 正确行为下的能量性质 |
|---|---|---|
| 0 | intersect 求交 | 产生信号：主光线求交，输出 G-buffer（命中标记/id/距离） |
| 1 | shade 着色 | Lambert 直接光照 + 全图能量重归一化，均值与上游一致 |
| 2 | accumulate 采样累积 | N_SAMPLES 次无偏抖动估计取均值，均值与上游一致 |
| 3 | denoise 降噪 | 3x3 钳位边界盒式滤波，权重归一化，均值与上游一致 |
| 4 | tonemap 色调映射 | 曝光归一化 + Reinhard，重标定回原均值尺度 |

每个阶段支持 `PASSTHROUGH`：输入 `memcpy` 到输出，不改一个像素。

## 用法

```
make
./pipeline <outdir> [cut]   # cut=0（默认）依次跑 cut1..cut5
```

`cut k` = 前 k 个阶段 ACTIVE，其余 PASSTHROUGH。每个 cut 输出：

- `cut{k}.pfm` — 线性 HDR 中间图（量测用）
- `cut{k}.ppm` — sRGB 8-bit 图（肉眼查看用）
- `cut{k}_luma.txt` — 各级之后的平均亮度

## 排查判据

从 cut1 到 cut5 逐级打开，比较相邻两级的全图平均亮度：
落差 > 1% 的那一级即为归一化写错的环节。
