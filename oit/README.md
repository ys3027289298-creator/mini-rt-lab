# oit

`make test` 跑自检。

- `oit_exact_sorted` — 按深度排序后精确 "over" 合成（线性算子的交换律形式）。
- `oit_weighted` — 深度指数加权累积，与输入顺序无关（逐位一致）。
- `oit_exact_weights` — 每层精确前景权重 `f_i = a_i * prod_{j<i}(1-a_j)` 与总透射率。
- `oit_kappa_for` — 给出使 `w_i <= f_i` 对全栈成立的最小 `kappa`。
- `oit_error_bound` — 加权结果与精确排序结果之差的严格上界。

零层（`n <= 0`）一律返回 `OIT_E_LAYERS`，不写输出。
