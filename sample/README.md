# sample

- `make test`：统一随机源自检（`RT_UNIFIED_RNG`，默认）。
- `make test-legacy`：关闭统一来源，走原始 Halton/分层路径。
- `make verify_hist`：构建直方图验收程序；运行
  `./verify_hist <outdir>` 生成 `histogram.csv`（stream 1）。
