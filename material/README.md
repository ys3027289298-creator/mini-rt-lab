# material

`make test` 跑自检。

含 Lambert / GGX 求值与余弦加权采样，以及次表面散射的扩散近似：

- 按扩散剖面（偶极子，匹配边界）做径向积分，`sss_capture()` 用闭式解给出半径 R 内的出射能量占比
- 半径由平均自由程定：剖面只依赖 `R/mfp`，`sss_radius_for_capture()` 反解给定占比对应的半径
- 能量损失：`sss_reflectance()` 把入射能量拆成侧向分辨项、Lambert 项与吸收项，三项恒等于 1
- 半径为零时严格退化为 Lambert 漫反射
- `sss_init()` 遇到非法参数（albedo 不在 [0,1]、mfp 非有限或 ≤ 0）返回出错的参数名
- 与同场景蒙特卡洛参考比对，误差上界见自检里的 `sss_mc_error_within_bound`
