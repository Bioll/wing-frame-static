"""读 Abaqus 四面体 .inp, 提取外表面渲染。确认哪些体被填充。"""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from collections import defaultdict

inp = sys.argv[1] if len(sys.argv) > 1 else "airfoil_frame_tet.inp"
nodes = {}; tets = []; sec = None; etype = "?"
for line in open(inp, encoding="utf-8", errors="ignore"):
    s = line.strip()
    if not s: continue
    if s.startswith("*"):
        u = s.upper()
        sec = "N" if u.startswith("*NODE") else ("E" if u.startswith("*ELEMENT") else None)
        if sec == "E":
            for tok in s.split(","):
                if "TYPE=" in tok.upper(): etype = tok.split("=")[1].strip()
        continue
    p = [x for x in s.replace(",", " ").split() if x]
    if sec == "N":
        nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
    elif sec == "E":
        c = [int(x) for x in p[1:]]
        if len(c) >= 4: tets.append(c[:4])
print(f"type={etype} nodes={len(nodes)} tets={len(tets)}")

# 外表面 = 只出现一次的三角面
fc = defaultdict(int); rep = {}
for t in tets:
    for tri in ((t[0],t[1],t[2]),(t[0],t[1],t[3]),(t[0],t[2],t[3]),(t[1],t[2],t[3])):
        k = tuple(sorted(tri)); fc[k]+=1; rep[k]=tri
bnd = [rep[k] for k,v in fc.items() if v==1]
print("boundary tri faces:", len(bnd))

ids = sorted(nodes); idx={n:i for i,n in enumerate(ids)}
xyz = np.array([nodes[n] for n in ids])
verts = [[xyz[idx[n]] for n in f] for f in bnd]
fig = plt.figure(figsize=(15,8))
for k,(e,a) in enumerate([(20,-60),(20,120)]):
    ax = fig.add_subplot(1,2,k+1,projection="3d")
    ax.add_collection3d(Poly3DCollection(verts, facecolor="#7fb27f",
                        edgecolor=(0,0,0,0.2), linewidths=0.1))
    mn,mx = xyz.min(0),xyz.max(0)
    ax.set_xlim(mn[0],mx[0]);ax.set_ylim(mn[1],mx[1]);ax.set_zlim(mn[2],mx[2])
    ax.set_box_aspect(mx-mn); ax.view_init(e,a); ax.set_title(f"view {k+1}")
fig.suptitle(f"airfoil frame TET mesh ({etype}, {len(tets)} tets)")
plt.tight_layout(); plt.savefig("airfoil_frame_tet.png", dpi=110)
print("wrote airfoil_frame_tet.png")
