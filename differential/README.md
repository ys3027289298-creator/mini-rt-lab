# differential

反射光线微分传播（Igehy 1999），给出纹理过滤足迹。`make test` 跑自检。

## 接口

| 函数 | 作用 |
|---|---|
| `rd_transfer` | 光线微分传播到命中曲面，切平面投影，足迹超上限时截断并返回 `RD_WARN_CLAMP` |
| `rd_reflect` | 镜面反射方向及其微分；曲率经 `dN = ±dP/R` 进入反射导数，法线朝向不影响结果 |
| `rd_measure` | 标量足迹宽度 `w = 0.5(|dP/dx|+|dP/dy|)` 与两轴长度 |
| `rd_uv_widths` | 解 2×2 Gram 系统求纹理空间过滤宽度（掠射角下分量除法会失效，这里不会） |

状态码：`RD_OK / RD_WARN_CLAMP / RD_ERR_NAN / RD_ERR_NEG / RD_ERR_PARAM / RD_ERR_DEGEN`。

## 关键公式

```
transfer:  dt = -((dP + t·dD)·N) / (D·N)
           dP' = dP + t·dD + dt·D
reflect:   r  = D − 2(D·N)N
           dr = −2[(D·dN)N + (D·N)dN]        （transfer 后入射 dD = 0）
```

注意 transfer 分子必须用曲面约束形式 `(dP+t·dD)·N`；Igehy 原文的
`dP·D + t·dD·D` 形式只在入射微分已切向时等价，pinhole 光线（dO=0）
下会静默跳过切平面投影——这正是朴素 `|dD|·t` 缩放在掠射角低估足迹的原因。

## 自检覆盖

- 平面镜：曲率项为零，二次传播后足迹宽度精确保持 1
- 单位球球心：解析宽度精确为 3（曲率使足迹放大 3 倍）
- 离轴球面：与 double 精度逐像素中心差分比对，误差上界
  `4/3|C−F| + 浮点界`，断言实测误差落在界内
- pinhole 地面：朴素距离缩放低估 ~5.3%，本模块给出精确值 3
- 极端曲率：返回 `RD_WARN_CLAMP`、宽度截断到 cap、光线方向不动；cap≤0 拒绝
- NaN/Inf、负平方长度、t<0、半径≤0、退化 UV 基全部拒绝
