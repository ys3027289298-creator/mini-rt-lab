# hybrid

烘焙光照 + 实时直接光的混合，C 实现，`make test` 跑自检并输出指标。

下文「纯实时」指逐像素精确阴影的直接光（实时渲染的质量基准），
不是无阴影直出 —— 与无阴影直出比会和被检查的缺陷同源。

## 问题

朴素叠加：`烘焙直接光 + 烘焙间接光 + 实时直接光`。

- 烘焙直接光被重复计入；
- 实时直接光没有遮挡判定，阴影里出现漏光。

## 混合

每个像素 `p`、每盏灯 `l`：

```
D_l(p) = I_l * max(0, N·L) / (|L-p|^2 + 1)     实时直接光（无遮挡）
V_l(p) = bilinear(烘焙可见性)                   烘焙可视性
λ_l(p) = 4 V (1 - V)                            漏光比例（插值不确定度）
Vc_l   = clamp(V * (1 - C_LEAK * λ), 0, 1)      按漏光比例补偿
L_hyb  = bake_indirect(p) + Σ_l D_l(p) * Vc_l(p)
```

- 不再叠加烘焙直接光，只保留烘焙间接光；
- 直接用烘焙可视性修正；
- `λ` 在 0/1 处为 0、在 0.5 处最大，正好对应双线性插值跨阴影边界的漏光带；补偿量与 `λ` 成正比；
- 烘焙数据缺失（双线性四个环绕纹素任一 `valid=0`）：该像素退回纯实时（逐像素精确直接光 + 逐像素精确间接光），并在 `scene.fallback` 标注。

## 指标

| 名称 | 含义 |
|---|---|
| `naive_leak_pixels` / `hybrid_leak_pixels` | 参考渲染中完全阴影（直接光≈0）但输出直接光 > `LEAK_TAU` 的像素数 |
| `leak_reduction_pct` | 相对朴素叠加的漏光像素下降百分比 |
| `leak_energy_reduction_pct` | 漏光能量（超出 `LEAK_TAU` 的部分求和）下降百分比 |
| `*_over_pixels` | 反向误差：参考中受光、但输出直接光比参考低 > `OVER_TAU` 的像素 |
| `*_max_err` / `*_mean_err` | 与纯实时参考渲染的误差 |
| `hybrid_err_max_vs_realtime` | 混合结果对纯实时参考的最大逐像素误差 |
| `hybrid_error_bound_max` | 逐像素三角不等式上界 |

上界形式：

```
valid:    |B_i - I| + Σ_l |D_l Vc_l - D_l s_l|
fallback: 0   (与参考逐像素一致)
```

其中 `B_i` 是采样到的烘焙间接光，`I` 是参考间接光，`s_l` 是参考可见性（0/1）。
自检断言实测误差不超过该上界（`error_bound_holds`），漏光像素至少减半，
烘焙缺失用例断言回退像素与纯实时参考误差为 0、标注掩码计数一致、非回退像素不受影响。

## 实测（默认场景，64×64）

```
shadow_pixels=920
naive:  leak_pixels=920  leak_energy=442.92  mean_err=0.8714  max_err=2.9507
hybrid: leak_pixels=252  leak_energy= 14.00  mean_err=0.0892  max_err=1.3644
        leak_reduction_pct=72.61  leak_energy_reduction_pct=96.84
        hybrid_err_max_vs_realtime=1.3644  hybrid_error_bound_max=1.3644
fallback（2 个纹素失效）：fallback_pixels=512  与纯实时误差=0.000000
```

## 补偿强度权衡

`make sweep` 扫描 `C_LEAK ∈ [0,1]`（`g_c_leak`）：

```
c_leak,leak_pixels,over_pixels,leak_energy,mean_err,max_err
0.0,554,833,86.8033,0.0604,1.2306
0.5,485,1045,50.4337,0.0735,0.9671
1.0,252,1138,14.0005,0.0892,1.3644
```

可见性遮罩本身就把漏光从 920 降到 554；补偿再降到 252，但会把受光侧的阴影边界压暗
（`over_pixels` 上升）。默认 `C_LEAK=1`，需要保守时调低即可。

## 调试图

`make images` 输出 `out/{reference,naive,hybrid,fallback_mask}.ppm`（PPM，可直接转 PNG），
回退像素在掩码中为品红。
