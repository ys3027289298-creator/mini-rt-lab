# mini-rt-lab

软件光线追踪基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `bvh/` | 包围盒树构建、中位数分裂、最近命中 |
| `sample/` | 分层与低差异采样、方差估计 |
| `material/` | Lambert 与 GGX 求值、余弦加权采样 |
| `tests/` | 自建最小测试入口（无第三方库）与变异缺陷注入，见 `tests/README.md` |

## 测试

```sh
make -C tests test            # 几何求交 + 材质求值两条路径的失败用例
make -C tests check-baseline  # 三个 demo 成功路径输出逐字节一致
make -C tests mutation        # 注入 24 处缺陷，捕获率 23/24 = 95.8%
```
