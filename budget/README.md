# budget

按块统计误差、导出热图、按分位数排名最差块、按预算分配样本。

`make test` 跑自检，并在当前目录写出 `error_heatmap_uniform.ppm`
（均分预算）和 `error_heatmap_opt.ppm`（按 √var 分配）。

自检要点：

- 热图峰值位置就是数值最大块的索引：导出后回读 PPM 找最亮像素，
  再和 `field_worst()` 对比，布局或打包写错会直接失败。
- 均分与最优分配用同一总预算，打印 `opt_over_uniform` 比值。
- 预算耗尽后仍超阈值的块按误差从大到小列出。
