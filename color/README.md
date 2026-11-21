# color

线性光渲染管线：解码到线性、在线性空间累积、最后统一编码输出。

- `Transfer` 描述传输函数。非线性类型（`TF_GAMMA`、`TF_SRGB`）必须带齐参数，
  缺参数或参数非法返回 `E_PARAM`，不会回退到默认值。
- `Pipeline` 的 `pipeline_add` 在入口把编码样本转成线性再累加；
  `pipeline_output` 归一化后只编码一次。
- `pixelwise_reference` 是逐像素精确转换的参考实现，
  自检要求管线结果与它一致到 1e-12。
- `roundtrip_bound` 给出解析上界（编码/线性/ΔE76），
  `measure_roundtrip` 在分段接缝处加密采样测得实际误差，自检要求测量值不超过上界。

`make test` 跑自检。
