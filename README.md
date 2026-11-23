# mini-rt-lab

软件光线追踪基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `bvh/` | 包围盒树构建、中位数分裂、最近命中 |
| `sample/` | 分层与低差异采样、方差估计 |
| `material/` | Lambert 与 GGX 求值、余弦加权采样 |
| `kernel/` | 统一随机源：衰减/材质/采样三站点共用单一随机流，含改动前后验收证据 |

`make test` 会先跑三个模块自检，再跑 `kernel/` 的验收套件；验收数据见
`kernel/evidence/REPORT.md`。
