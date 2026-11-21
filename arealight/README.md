# arealight

`make test` 跑自检。

面光源按网格格心采样，替代"只取中心"：

- **格心采样**：把光源矩形切成 `nx*ny` 个格子，每个格心按解析软影权重 `phi` 取值：
  完全遮挡一端权重 0，完全可见一端权重 1，半影里按 smoothstep 过渡。网格是沿
  遮挡边坐标 `q = x - M*y` 的一维中点积分，`phi` 正是可见性权重。
- **可见性权重**：每个格心只算它"看到"的光——遮挡体在该 `q` 处把光切成前缀，
  可见比例就是 `phi((q-e)/W)`，`e` 是遮挡边投到光源平面的位置。光源与直线遮挡边
  都是 `q` 的等值线，所以积分降到一维；`rho(q)` 是光源在 `q` 处的弦长（体积因子）。
- **网格密度调噪声**：传入种子时每个样本在格子内按哈希抖动，把偏差变成噪声，
  标准差有界 `sigma <= 0.5*sqrt(4n+2)/n^2`。`arealight_grid_for_noise(eps)`
  返回满足目标 `eps` 的最小网格。
- **误差上界**：每个格子累加一个严格上界——跨越半影边的格子按 `phi` 在格内的
  变化量计，其余格子按中点求积误差 `(1/2)*Var(q)*|phi''|` 计，
  `Var(q) = (hx^2 + M^2 hy^2)/12`。自检验证该上界始终包住实际误差，且随网格加密
  单调下降。
- **零尺寸网格**：任一维为 0 时返回 `SHADOW_NO_GRID`，明确报告"无法切分"，
  不再退化成单点采样。

## 接口

```c
ShadowStatus arealight_shadow(V3 P, int nx, int ny,
                              const double* jitter_seed, ShadowResult* out);
int arealight_grid_for_noise(double eps);
```

`ShadowResult` 含 `vis`、`grad`、`samples`、`cell_q`、`err_v`、`err_g`、
`noise_sigma`。`jitter_seed == NULL` 时按格心确定采样，`noise_sigma = 0`。

## 自检

```
$ make test
ref_matches_dense_quadrature=1
ref_gradient_matches_fd=1
refine_lowers_error_monotone=1
bound_monotone=1
bound_contains_error=1
zero_grid_reported=1
behind_occluder_reported=1
center_only_error_large=1
noise_target_attained=1
noise_bound_monotone=1
noise_bound_holds_empirically=1
jitter_mean_tracks_reference=1
all_ok=1
```

注：上界与网格加密单调对应的是**认证上界**；个别采样点在换成恰好不利的网格相位时
实际误差会轻微回升（走样），自检里用整片接收点的最大误差来验证单调下降。
