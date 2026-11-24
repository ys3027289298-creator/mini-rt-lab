# tools

五阶段管线直通排查的离线分析工具，仅依赖 Python 标准库，用 `python3 -I` 运行。

## analyze_pipeline.py

直接解析 `pipeline` 写出的 PFM 中间图（不读程序自报的亮度文本，量测独立），
重算 cut1..cut5 的全图平均亮度，按 1% 阈值标记落差超标的阶段，并生成
SVG 亮度曲线 + HTML 报告（含每个阶段的落差数据表和五段中间图）。

```
python3 -I analyze_pipeline.py artifacts/runs/buggy artifacts/runs/fixed \
        --out artifacts/pipeline_report.html
```

## render_check.py

无浏览器环境下的目视检查：自带 5x7 点阵字体和 zlib PNG 编码器，
把 PPM 中间图转 PNG，并把亮度曲线栅格化为 `brightness_curve.png`。

```
python3 -I render_check.py artifacts/runs/buggy artifacts/runs/fixed
```

## 复现流程

```
cd pipeline
make clean && make all pipeline-buggy
./pipeline-buggy ../artifacts/runs/buggy     # 缺陷版：累积/降噪各掉 50%
./pipeline         ../artifacts/runs/fixed   # 修复版：曲线平坦
cd ..
python3 -I tools/analyze_pipeline.py artifacts/runs/buggy artifacts/runs/fixed
python3 -I tools/render_check.py     artifacts/runs/buggy artifacts/runs/fixed
```
