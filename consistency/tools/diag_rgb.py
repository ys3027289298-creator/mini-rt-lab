import sys
import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from consist_check import read_fbuf

r = read_fbuf(sys.argv[1])
t = read_fbuf(sys.argv[2])
rid = r[:, :, 4].round().astype(int)
tid = t[:, :, 4].round().astype(int)
for tri in range(3):
    m = (rid == tri) & (tid == tri)
    if not m.any():
        continue
    rc = r[m][:, 0:3]
    tc = t[m][:, 0:3]
    print(f"--- tri {tri} ({m.sum()} px) ---")
    print("  rast unique-ish rgb mean/min/max:")
    print("   ", rc.mean(0).round(5), rc.min(0).round(5), rc.max(0).round(5))
    print("  trace rgb mean/min/max:")
    print("   ", tc.mean(0).round(5), tc.min(0).round(5), tc.max(0).round(5))
    # expected flat: albedo * 0.15
