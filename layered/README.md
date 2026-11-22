# layered

多层平面介质栈（多片玻璃叠加）的颜色/能量传输求解，单文件 C，`make test` 跑自检。

## 模型

介质按序排列 `media[0] | media[1] … | media[n-1]`，首尾为半无限空间，中间介质各带
物理厚度。每个 RGB 通道独立使用该通道的折射率 `eta` 与吸收系数 `sigma`（波长相关 IOR
与吸收正是修正"多层玻璃偏色"的关键——不能把每层 albedo 直接相乘）。

- **界面**：Snell 定律折射/反射；非偏振 Fresnel 强度系数（s/p 平均）分配能量，
  正反方向互易；全反射（TIR）时反射率为 1、透射为 0。
- **层间衰减**：Beer–Lambert，路径长度按几何距离 `d / cos_theta`，
  `tau = exp(-sigma * d / cos_theta)`。

### 两个求解器

1. `layered_transport`：**解析非相干多层透射率**。净辐射（net-radiation）递推，从最
   后一个界面反向合并，自动包含所有层间多次反射，作为参考基准。厚玻璃表面间距远大于
   相干长度，非相干强度叠加是正确模型。
2. `layered_trace`：**确定性路径分裂追踪**。在每个界面按 Fresnel 分裂为折射/反射两条
   子路径，深度超过 `max_depth` 的路径权重计入 `unresolved`。

### 保证（由自检覆盖）

- **误差上界**：截断只丢路径不造能量，故
  `|T_trace - T_analytic| <= unresolved`，且为单侧界（截断只减少透射）；
  无损单层上该上界随深度按几何级数（公比 F²）收敛到零。
- **能量不超过入射**：吸收在每次衰减发生时即时记账，恒有
  `T + R + A + unresolved = 1`（追踪器）与 `T + R <= 1`（解析器）。
- **零层报告**：介质数 < 2 或 `nslab == 0` 返回 `LAYERED_EMPTY`，不静默当作真空；
  非法 IOR/厚度/入射角返回 `LAYERED_PARAM`，超容量返回 `LAYERED_CAP`。

## API 摘要

```c
int layered_from_media(Stack*, const Medium* media, const double* thick, int nmedia);
int layered_from_slabs(Stack*, const Medium* ambient, const Medium* gap, double gap_d,
                       const Medium* slabs, const double* slab_d, int nslab);
int layered_transport(const Stack*, double cos_incident, Transport* out); /* T,R,A 解析 */
int layered_trace(const Stack*, double cos_incident, int max_depth, TraceResult* out);
```
