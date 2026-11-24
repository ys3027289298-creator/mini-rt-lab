# material

`make test` 跑自检（Lambert/GGX 求值、余弦采样、VNDF 采样恒等式、
掠射 Taylor 模型、α=0 镜像退化）。

`make exp` 跑掠射高光实验（`src/highlight.c`）：平面 + 近地平线高斯
环境光，对比余弦采样 vs VNDF + Taylor 掠射带替换，多次重复测方差，
产物（PFM 帧、裁图 PNG、对比 SVG、stats.json）写入 `../artifacts/`，
分析见 `../artifacts/highlight_report.md`。

- `src/ggx.h` — GGX 全家桶（D/G1/G2/Schlick）、Heitz VNDF 采样、
  q(t)=pdf·brdf 在掠射方向的二阶 Taylor 展开，供自检与实验共用。
- `tools/compose.py` — 裁图合成与方差标注（纯标准库）。
