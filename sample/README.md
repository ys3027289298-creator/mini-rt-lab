# sample

`make test` 跑自检。

`cmj_layer(out, grid, block, layer)`：相关多抖动采样。一层是 grid×grid 层叠拼块、每格一个样本；连续样本按 `block` 分块，块内共享同一抖动向量（层内相关），块间各自独立抽抖动，不同 `layer` 重新播种（层间独立）。任意块尺寸下样本都不越出自己的格子，因此二维直方图均匀性有确定性上界（`hist2d_bound`，k 整除 grid 时误差恒为 0）。`block < 1` 或 `block > grid*grid` 返回 -1 拒绝。自检打印块尺寸-相关性数据并与纯随机序列比对直方图误差。
