import subprocess
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from consist_check import read_fbuf

SCENES = ["clean", "lab"]
SIZES = [80, 160, 320, 640]


def render(impl, scene, w, h, out):
    subprocess.run([f"./{impl}", f"scenes/{scene}.scene", str(w), str(h), out],
                   check=True, stdout=subprocess.DEVNULL)


for scene in SCENES:
    print(f"=== {scene} ===")
    for w in SIZES:
        h = int(w * 0.75)
        render("raster", scene, w, h, "out/_r.fbuf")
        render("tracer", scene, w, h, "out/_t.fbuf")
        r = read_fbuf("out/_r.fbuf")
        t = read_fbuf("out/_t.fbuf")
        rid = r[:, :, 4].round().astype(int)
        tid = t[:, :, 4].round().astype(int)
        same = (rid == tid) & (rid >= 0)
        idd = (rid != tid) & ((rid >= 0) | (tid >= 0))
        cd = np.abs(r[:, :, 0:3].astype(float) - t[:, :, 0:3].astype(float)).max(2)
        zr, zt = r[:, :, 3].astype(float), t[:, :, 3].astype(float)
        zs = np.maximum(np.abs(zr), np.abs(zt))
        dd = np.abs(zr - zt) / np.where(zs > 1e-6, zs, 1)
        cds = cd[same]
        dds = dd[same]
        print(f" {w}x{h}: sameid={same.sum():6d} iddiff={idd.sum():5d} | "
              f"color q50/90/99/max={np.quantile(cds,.5):.4f}/{np.quantile(cds,.9):.4f}/"
              f"{np.quantile(cds,.99):.4f}/{cds.max():.4f} | "
              f"depthrel max={dds.max():.2e}")
