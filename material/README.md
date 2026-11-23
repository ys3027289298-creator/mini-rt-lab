# material

- `make test`：统一随机源自检（`RT_UNIFIED_RNG`，默认）。
- `make test-legacy`：关闭统一来源，走原始网格算术采样路径。
- `make verify_variance`：构建方差验收程序；运行
  `./verify_variance <outdir>` 生成 `variance.json`（stream 2）。
