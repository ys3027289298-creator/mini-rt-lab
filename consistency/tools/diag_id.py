import subprocess
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from consist_check import read_fbuf

w, h = 320, 240
subprocess.run(["./raster", "scenes/lab.scene", str(w), str(h), "out/_r.fbuf"],
               check=True, stdout=subprocess.DEVNULL)
subprocess.run(["./tracer", "scenes/lab.scene", str(w), str(h), "out/_t.fbuf"],
               check=True, stdout=subprocess.DEVNULL)
r = read_fbuf("out/_r.fbuf")
t = read_fbuf("out/_t.fbuf")
rid = r[:, :, 4].round().astype(int)
tid = t[:, :, 4].round().astype(int)
idd = (rid != tid) & ((rid >= 0) | (tid >= 0))
print("id-diff pixels:", idd.sum())
ys, xs = np.nonzero(idd)
for i in range(min(20, len(ys))):
    y, x = ys[i], xs[i]
    print(f"({x},{y}) rast_id={rid[y,x]} trace_id={tid[y,x]} "
          f"zr={r[y,x,3]:.5f} zt={t[y,x,3]:.5f}")
# per-id coverage counts to confirm all triangles actually rasterize
for name, ids in (("rast", rid), ("trace", tid)):
    u, c = np.unique(ids[ids >= 0], return_counts=True)
    print(name, dict(zip(u.tolist(), c.tolist())))
