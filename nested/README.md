# nested

玻璃球（半径 1，IOR 1.5）内含一颗偏移的 Lambertian 小球（半径 0.55）。
对比两种追踪器：

- `buggy`：复现遗留行为——外层介质表面折射后只再做**一次**求交；若该
  条光线穿出玻璃壳就直接终止。经全内反射（TIR）折返再向内的路径全部
  丢失，内层球边缘因此被切成块状缺口。
- `walk`：沿光路连续求交。入口、出口、TIR 反弹、内层 Lambertian 反弹
  都是同一条路径的下一段，直到吞吐量低于能量阈值或穿透次数达到上限。

## 终止参数如何从内层平均反射率 rho 推出（非经验值）

每穿透一次内层球，吞吐量乘一次 rho，K 次穿透后残余能量为 rho^K。
8 位输出无法区分低于一个量化步长 LSB = 1/255 的贡献，因此：

```
K   = min { k >= 1 : rho^k <= LSB }      = ceil( ln(LSB) / ln(rho) )
tau = rho^K                               （此时恰有 tau <= LSB）
```

截断误差被严格限制在输出不可表示的能量以下。`nested.c` 的
`derive_limits()` 实现该推导，自检校验了两个反照率下的整数解：

| rho | K（穿透上限） | tau = rho^K（能量阈值） |
|---|---|---|
| 0.25 | 4 | 0.003906 |
| 0.75 | 20 | 0.003171 |

反照率越高、能量衰减越慢，需要的穿透次数越多——上限由 rho 决定，不靠调参。

## 用法

```
make test                       # 几何/Fresnel/TIR/推导/RNG 自检
./nested W H SPP rho buggy|walk K tau seed out.ppm
                                # K=-1 / tau=-1 表示按 rho 自动推导
```

每次渲染同时输出颜色图 `out.ppm` 和内层命中掩码 `out.mask.ppm`。

## 文件

| 路径 | 内容 |
|---|---|
| `src/nested.c` | 两种追踪器、Fresnel/Schlick 之外的精确反射率、TIR、rho 推导、自检 |
| `tools/analyze.py` | 以 walk 推导参数渲染为参考，统计边缘带缺口/差异块与亮度偏置 |
| `run_experiments.sh` | 两个反照率的完整实验（见 RESULTS.md） |
| `renders/exp/` | PPM/PNG 结果、边缘裁剪、stats.json/stats.txt、run.log |
