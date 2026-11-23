# bvh

- `make test`：统一随机源自检（`RT_UNIFIED_RNG`，默认），构建前对图元
  顺序做一次来自统一流的 Fisher-Yates 置换。
- `make test-legacy`：关闭统一来源，走原始中位数分裂路径。
- `make verify_image`：构建图像验收程序；运行
  `./verify_image <outdir>` 生成 `image.ppm` 与 `image_hash.json`
  （stream 3，FNV-1a64 哈希）。
