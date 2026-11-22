# atlas

贴图图集打包：把多张小贴图按面积降序做 shelf 装箱进一张大图集，
渲染时一次绑定代替逐张贴图绑定，降低绘制批次。

## 构建与自测

```
make test
```

## 自测内容

- `pack_ok` / `all_placed` / `no_overlap`：全部贴图装箱成功且无重叠
- `zero_size_rejected`：宽或高为 0 的贴图拒绝装箱（返回错误，不参与布局）
- `pixel_identical`：同一场景分别用逐张贴图与图集渲染，帧缓冲逐字节一致
- `batches_reduced`：图集批次 < 逐张贴图批次
- `utilization_above_half`：图集面积利用率 > 50%

## 数据输出

`make test` 打印 `utilization` / `waste` / `batches_per_texture` /
`batches_atlas` / `batch_reduction` 等 key=value 数据行。

## 画面对照

运行后生成 PPM（可用任意看图工具或 `convert x.ppm x.png` 查看）：

- `atlas.ppm`：打包后的图集（黑色为 gutter 与未利用区域）
- `render_individual.ppm` / `render_atlas.ppm`：两种方式的渲染结果
- `compare.ppm`：左逐张贴图、右图集并排对照，底部为 ×16 放大的差异条
  （逐像素一致时为纯黑）

## 说明

- 采样为最近邻，texel 索引两条路径完全一致，因此渲染结果逐像素相同；
  若换双线性过滤，需要 gutter 内扩边（duplicate edge texels）才能避免串色。
- 装箱宽度在 `[sqrt(总面积), 2*sqrt(总面积)]` 内扫描，取图集面积最小者。
- shelf 装箱简单确定，但利用率一般（本例 ~55%）；要更高利用率可换
  MaxRects 等启发式，接口不变。
