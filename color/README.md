# color

贴图色彩解码一致性自检，`make test` 跑自检。

管线约定：纹素先按 sRGB EOTF 解码到线性空间，过滤/混合/着色都在线性空间做，
仅在输出显示时编码回非线性空间。

| 检查 | 含义 |
|---|---|
| `decode_*` | sRGB 分段解码（IEC 61966-2-1）端点、已知值、分段点连续性、单调性 |
| `decode_vs_double_ref` | float 解析解码与 double 参考实现的最大偏差 |
| `lut_pixel_exact` | 8-bit LUT 采样路径与解析解码 256 级逐像素一致 |
| `srgb_rt_*` | 非线性→线性→非线性往返误差，float 与 8-bit 分别给实测最大值和上界 |
| `blend_matches_analytic` | 线性空间双线性混合与解析路径逐点一致 |
| `gammaspace_*` | 对照实验：sRGB 空间混合后解码只会更暗（EOTF 凸性），量化贴图发暗原因 |
| `gamma_*_bound` | 通用 gamma(2.2) 往返误差上界 |
| `gamma_invalid_reports_value` | gamma 取 0/负/NaN/Inf 时拒绝并报出具体取值 |
