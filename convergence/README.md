# convergence

两种低差异序列 (Halton vs Sobol) 的收敛指数对比。

`make test` 跑自检, `make report` 输出对比报告, `./converge k` 指定样本数档位个数
(默认 9 档, N = 16..4096; k < 4 时各序列报告"不可拟合")。

方法: 每个 N 用 Cranley-Patterson 随机移位重复 32 次取 RMSE,
在 log-log 下最小二乘拟合 `err = C * N^-alpha`, 报告:

- 指数 alpha 及 95% 置信区间 (t 分布);
- 同样本 (最大 N) 下两序列的误差比例;
- 残差上界: 所有观测误差 <= 拟合值 x exp(max|log 残差|);
- 样本点不足 (< 4)、误差非正或样本数无变化时报"不可拟合"及原因。
