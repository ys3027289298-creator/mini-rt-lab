# mini-rt-lab

软件光线追踪基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `bvh/` | 包围盒树构建、中位数分裂、最近命中 |
| `sample/` | 分层与低差异采样、方差估计 |
| `material/` | Lambert 与 GGX 求值、余弦加权采样 |

## 无依赖构建（每个子目录独立）

只依赖 C11 编译器和 libm，无第三方库。三个目录互不依赖，可分别构建：

```sh
make -C bvh            # 或: cd bvh && make
make -C sample
make -C material
make -C bvh test       # 自检；render 目标生成演示图 bvh.ppm
```

默认编译选项：`-O2 -std=c11 -Wall -Wextra -Werror -ffp-contract=off`（零告警，告警即构建失败）。

## 渲染

每个程序支持统一的渲染子命令：

```sh
./bvh render out.ppm [--size N] [--strict|--no-strict]
```

- `--size N`：图像边长，整数 1..1024，默认 128；非法参数（非数字、越界、缺值、未知选项、输出路径不可写）打印错误说明并以非零码退出。
- `--strict`（默认）/ `--no-strict`：开关渲染期运行时校验；只校验不改值，两种模式产物逐字节相同。

## 可复现性

渲染路径只用整数运算、IEEE 浮点四则运算和 `sqrtf`（IEEE 754 规定正确舍入），不使用 `cosf`/`sinf` 等 libm 超越函数（其结果因 libm 实现而异），并显式 `-ffp-contract=off` 禁止 FMA 收缩。因此同一源码在不同编译器、优化级别、机器上渲染结果逐像素一致。

一键复现全部证据（写入 `artifacts/repro/`，任一项失败即中止）：

```sh
./repro.sh
```

覆盖：两种编译选项（gcc -O2 / clang -O3）零告警构建、自检输出一致、渲染结果逐像素比较（`tools/ppmdiff.c`，diff 必须为 0）、8 类非法参数报错、strict 开关产物逐字节一致。
