# sample

无依赖构建（仅需 C11 编译器 + libm）：

```sh
make            # 构建，-Wall -Wextra -Werror 零告警
make test       # 自检
make render     # 生成 sample.ppm（./sample render sample.ppm [--size N] [--strict|--no-strict]）
```

跨编译选项/跨机器的逐像素可复现验证见仓库根目录 `repro.sh`。
