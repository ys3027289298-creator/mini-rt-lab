# mini-rt-lab

软件光线追踪基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `bvh/` | 包围盒树构建、中位数分裂、最近命中 |
| `sample/` | 分层与低差异采样、方差估计 |
| `material/` | Lambert 与 GGX 求值、余弦加权采样 |
| `perf/` | 三处渲染路径的统一性能预算、低/高两档验收（耗时表、峰值内存、图像哈希） |

## 性能预算验收（两档）

`perf/` 用一份 `perf/budget.conf` 管住 `bvh/`、`sample/`、`material/`
三条路径，低/高两档都必须给出耗时表、峰值内存、图像哈希与复现命令，
且都不得超过上界；关闭预算检查后输出逐像素不变；任一一档数据缺失即判未完成。

```sh
make -C perf accept     # 自检 + 两档 + 校验 + 关闭检查逐像素不变
```

详见 [`perf/README.md`](perf/README.md)。产物在 `artifacts/perf/`
（`REPORT.md`、每档 `metrics.json` 与三张 PPM）。
