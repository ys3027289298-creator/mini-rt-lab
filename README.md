# mini-rt-lab

软件光线追踪基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `bvh/` | 包围盒树构建、中位数分裂、最近命中 |
| `sample/` | 分层与低差异采样、方差估计 |
| `material/` | Lambert 与 GGX 求值、余弦加权采样 |
| `convergence/` | Halton vs Sobol 收敛指数拟合对比（置信区间、误差比、残差上界） |
