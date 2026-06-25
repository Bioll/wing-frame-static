"""读 Abaqus .inp (节点+tri/quad单元) 渲染成 png。用于看 HyperMesh 网格。"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

inp = sys.argv[1] if len(sys.argv) > 1 else "airfoil_hm.inp"
png = sys.argv[2] if len(sys.argv) > 2 else "airfoil_hm.png"

nodes = {}
faces = []      # 每个单元的节点id列表(3或4)
sec = None
for line in open(inp, encoding="utf-8", errors="ignore"):
    s = line.strip()
    if not s:
        continue
    if s.startswith("*"):
        u = s.upper()
        if u.startswith("*NODE"):
            sec = "N"
        elif u.startswith("*ELEMENT"):
            sec = "E"
        else:
            sec = None
        continue
    p = [x for x in s.replace(",", " ").split() if x]
    if sec == "N":
        nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
    elif sec == "E":
        conn = [int(x) for x in p[1:]]
        if len(conn) in (3, 4):
            faces.append(conn)

ids = sorted(nodes)
idx = {n: i for i, n in enumerate(ids)}
xyz = np.array([nodes[n] for n in ids])
print(f"nodes={len(ids)} elems={len(faces)}")

verts = [[xyz[idx[n]] for n in f] for f in faces]
fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111, projection="3d")
pc = Poly3DCollection(verts, facecolor=(0.9, 0.78, 0.55),
                      edgecolor=(0.15, 0.15, 0.15), linewidths=0.15)
ax.add_collection3d(pc)
mn, mx = xyz.min(0), xyz.max(0)
ax.set_xlim(mn[0], mx[0]); ax.set_ylim(mn[1], mx[1]); ax.set_zlim(mn[2], mx[2])
ax.set_box_aspect(mx - mn)
ax.view_init(elev=20, azim=-60)
ax.set_title(f"airfoil mesh (HyperMesh, {len(faces)} elems)")
plt.tight_layout()
plt.savefig(png, dpi=110)
print("wrote", png)
