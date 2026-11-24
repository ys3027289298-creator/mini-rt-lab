#!/usr/bin/env python3
"""analyze_pipeline.py — 直通排查的离线分析器（仅用 Python 标准库）。

输入：pipeline 程序写出的 run 目录，内含 cut1..cut5 的 PFM 中间图。
做三件事：
  1. 直接解析 PFM（不依赖管线自报的 luma 文本），独立重算每张中间图的
     全图平均亮度，得到"逐级打开"曲线；
  2. 计算相邻阶段的亮度落差，> 1% 即标记为归一化嫌疑环节；
  3. 输出 SVG 亮度曲线 + HTML 报告（含数据表，对比度不足的颜色用标签兜底）。

用法：
  python3 -I analyze_pipeline.py <run_dir> [<run_dir> ...] [--out report.html]
"""
import argparse
import html
import os
import struct
import sys

STAGES = ["intersect", "shade", "accumulate", "denoise", "tonemap"]
STAGE_CN = ["求交", "着色", "采样累积", "降噪", "色调映射"]
DROP_THRESHOLD = 0.01

# dataviz 技能默认分类色板（已通过 validate_palette.js，light 模式）
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
DANGER = "#e34948"


def read_pfm(path):
    """返回 (width, height, [(r,g,b), ...])，小端 PF 三通道 float。"""
    with open(path, "rb") as f:
        magic = f.readline().strip()
        if magic != b"PF":
            raise ValueError(f"{path}: 不是 PFM RGB 文件")
        # 头部允许注释行
        dim = f.readline()
        while dim.startswith(b"#"):
            dim = f.readline()
        w, h = map(int, dim.split())
        scale_line = f.readline().strip()
        little = scale_line.startswith(b"-")
        count = w * h * 3
        data = f.read(count * 4)
        if len(data) != count * 4:
            raise ValueError(f"{path}: 数据截断")
        fmt = "<" if little else ">"
        vals = struct.unpack(f"{fmt}{count}f", data)
        px = [(vals[i], vals[i + 1], vals[i + 2]) for i in range(0, count, 3)]
        return w, h, px


def mean_luma(px):
    return sum((r + g + b) / 3.0 for r, g, b in px) / len(px)


def analyze_run(run_dir):
    """每个 cut 给一张最终 PFM；cut k 的终图 = 前 k 级 ACTIVE 后的画面。"""
    points = []  # [(stage_index_just_opened, cut, luma)]
    for cut in range(1, 6):
        path = os.path.join(run_dir, f"cut{cut}.pfm")
        w, h, px = read_pfm(path)
        points.append((cut - 1, cut, mean_luma(px), (w, h)))
    return points


def stage_deltas(points):
    """相邻 cut 之间的落差（在新打开的那一级上发生）。"""
    rows = []
    for idx in range(1, len(points)):
        prev = points[idx - 1][2]
        cur = points[idx][2]
        drop = (prev - cur) / prev if prev > 1e-12 else 0.0
        rows.append((points[idx][0], prev, cur, drop, drop > DROP_THRESHOLD))
    return rows


def build_svg(runs, width=860, height=420):
    """多条 cut 曲线：x=新打开阶段，y=平均亮度。"""
    ml, mr, mt, mb = 70, 24, 30, 64
    pw, ph = width - ml - mr, height - mt - mb
    all_y = [p[2] for r in runs for p in r["points"]]
    ymax = max(all_y) * 1.12
    ymin = 0.0

    def x_of(i):
        return ml + pw * i / 4.0

    def y_of(v):
        return mt + ph * (1 - (v - ymin) / (ymax - ymin))

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" font-family="system-ui,sans-serif">',
        f'<rect width="{width}" height="{height}" fill="{SURFACE}"/>',
        f'<text x="{width/2:.0f}" y="20" text-anchor="middle" font-size="15" '
        f'font-weight="600" fill="{INK}">逐级打开后的全图平均亮度（PFM 独立量测）</text>',
    ]
    # y 网格 + 刻度
    for g in range(5):
        v = ymin + (ymax - ymin) * g / 4.0
        yy = y_of(v)
        parts.append(f'<line x1="{ml}" y1="{yy:.1f}" x2="{ml+pw}" y2="{yy:.1f}" '
                     f'stroke="{GRID}" stroke-width="1"/>')
        parts.append(f'<text x="{ml-8}" y="{yy+4:.1f}" text-anchor="end" font-size="11" '
                     f'fill="{INK_2}">{v:.4f}</text>')
    # 轴
    parts.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="{AXIS}"/>')
    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="{AXIS}"/>')
    # x 标签：阶段名 + cut 编号
    for i, (cn, en) in enumerate(zip(STAGE_CN, STAGES)):
        xx = x_of(i)
        parts.append(f'<text x="{xx:.1f}" y="{mt+ph+22}" text-anchor="middle" '
                     f'font-size="12" font-weight="600" fill="{INK}">{cn}</text>')
        parts.append(f'<text x="{xx:.1f}" y="{mt+ph+38}" text-anchor="middle" '
                     f'font-size="10.5" fill="{INK_2}">{en} · cut{i+1}</text>')
    parts.append(f'<text x="18" y="{mt+ph/2:.0f}" transform="rotate(-90 18 {mt+ph/2:.0f})" '
                 f'text-anchor="middle" font-size="12" fill="{INK_2}">平均亮度 (linear luma)</text>')

    # 每条 run 一条折线
    for ri, run in enumerate(runs):
        color = SERIES_COLORS[ri % len(SERIES_COLORS)]
        coords = [(x_of(p[0]), y_of(p[2])) for p in run["points"]]
        d = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in coords)
        parts.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="2"/>')
        for (x, y), p in zip(coords, run["points"]):
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{color}" '
                         f'stroke="{SURFACE}" stroke-width="2"><title>{html.escape(run["name"])} '
                         f'cut{p[1]}: {p[2]:.6f}</title></circle>')
        # 直接标注系列名（对比度 WARN 的兜底：文字用墨色，不靠颜色辨认）
        lx, ly = coords[-1]
        parts.append(f'<text x="{lx-6:.1f}" y="{ly-8:.1f}" text-anchor="end" '
                     f'font-size="11" font-weight="600" fill="{INK}">{html.escape(run["name"])}</text>')

    # 1% 落差标记（仅第一条 run，即主分析对象）
    run = runs[0]
    for row in stage_deltas(run["points"]):
        if row[4]:
            stage_i, prev, cur, drop, _ = row
            x = x_of(stage_i)
            parts.append(f'<text x="{x:.1f}" y="{mt+14}" text-anchor="middle" font-size="11" '
                         f'font-weight="700" fill="{DANGER}">▼ {drop*100:.1f}%</text>')
    parts.append("</svg>")
    return "".join(parts)


def build_html(runs, out_path):
    svg = build_svg(runs)
    chunks = [
        "<!doctype html><html lang='zh'><head><meta charset='utf-8'>",
        "<title>五阶段直通排查报告</title>",
        "<style>body{font-family:system-ui,sans-serif;color:#0b0b0b;background:#fcfcfb;"
        "max-width:920px;margin:24px auto;padding:0 16px;line-height:1.5}"
        "table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13px}"
        "th,td{border:1px solid #c3c2b7;padding:6px 9px;text-align:right}"
        "th:first-child,td:first-child{text-align:left}"
        "th{background:#f1f0ea}.bad{color:#e34948;font-weight:700}"
        ".ok{color:#006300;font-weight:600}h2{margin-top:28px;font-size:17px}"
        "code{background:#f1f0ea;padding:1px 5px;border-radius:3px}</style></head><body>",
        "<h1>五阶段渲染管线 · 直通模式能量排查</h1>",
        svg,
        "<p>判据：相邻两级平均亮度落差超过 <b>1%</b> 即判定该级归一化有误；"
        "直通级必须零损耗（曲线上表现为与前一 cut 重合）。</p>",
    ]
    for run in runs:
        deltas = stage_deltas(run["points"])
        chunks.append(f"<h2>{html.escape(run['name'])} <small style='font-weight:400;color:#52514e'>"
                      f"({html.escape(run['dir'])})</small></h2>")
        chunks.append("<table><tr><th>新打开阶段</th><th>cut</th><th>上游亮度</th>"
                      "<th>本级亮度</th><th>落差</th><th>判定</th></tr>")
        base = run["points"][0]
        chunks.append(f"<tr><td>{STAGE_CN[0]}（基准）</td><td>cut{base[1]}</td>"
                      f"<td>—</td><td>{base[2]:.6f}</td><td>—</td><td class='ok'>基准</td></tr>")
        for stage_i, prev, cur, drop, bad in deltas:
            cls = "bad" if bad else "ok"
            verdict = "归一化嫌疑" if bad else "通过(&lt;1%)"
            chunks.append(
                f"<tr><td>{STAGE_CN[stage_i]} <code>{STAGES[stage_i]}</code></td>"
                f"<td>cut{stage_i+1}</td><td>{prev:.6f}</td><td>{cur:.6f}</td>"
                f"<td class='{cls}'>{drop*100:.2f}%</td><td class='{cls}'>{verdict}</td></tr>")
        chunks.append("</table>")
        total = (base[2] - run["points"][-1][2]) / base[2]
        chunks.append(f"<p>全链路累计变化：<b class='{'bad' if total > DROP_THRESHOLD else 'ok'}'>"
                      f"{total*100:.2f}%</b>（终值 {run['points'][-1][2]:.6f} / 基准 {base[2]:.6f}）</p>")
        # 五段直通中间图（PNG 由 render_check.py 生成，相对报告路径）
        rel = os.path.relpath(run["dir"], os.path.dirname(out_path))
        chunks.append("<div style='display:flex;gap:10px;flex-wrap:wrap'>")
        for cut in range(1, 6):
            src = f"{rel}/cut{cut}.png"
            chunks.append(
                f"<figure style='margin:0;flex:1;min-width:150px'>"
                f"<img src='{src}' width='100%' style='border:1px solid #c3c2b7' "
                f"alt='{html.escape(run['name'])} cut{cut} 中间图'>"
                f"<figcaption style='font-size:12px;color:#52514e;margin-top:4px'>"
                f"cut{cut} · {STAGE_CN[cut-1]}</figcaption></figure>")
        chunks.append("</div>")
    chunks.append("</body></html>")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("".join(chunks))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    runs = []
    for d in args.run_dirs:
        points = analyze_run(d)
        name = os.path.basename(os.path.normpath(d))
        runs.append({"dir": d, "name": name, "points": points})
        print(f"== {name} ({d}) ==")
        base = points[0][2]
        print(f"  cut1 基准 ({STAGES[0]}): {base:.6f}")
        for stage_i, prev, cur, drop, bad in stage_deltas(points):
            flag = "  <-- DROP > 1%" if bad else ""
            print(f"  cut{stage_i+1} {STAGES[stage_i]:<10}: {cur:.6f}  落差 {drop*100:6.2f}%{flag}")
        total = (base - points[-1][2]) / base
        print(f"  全链路累计: {total*100:.2f}%  终值/基准 = {points[-1][2]/base:.4f}")
    out = args.out or os.path.join(os.path.dirname(os.path.normpath(args.run_dirs[0])),
                                   "pipeline_report.html")
    build_html(runs, out)
    print(f"报告已写入 {out}")


if __name__ == "__main__":
    sys.exit(main())
