"""去掉 vol 5(碳棒)/6(舵机)/7(蒙皮), 保留 1-4 内部骨架。导出新STEP + 渲染确认。"""
import gmsh, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Patch

REMOVE = [5, 6, 7]
KEEP   = [1, 2, 3, 4]
OUT_STEP = "airfoil_frame.stp"
OUT_PNG  = "airfoil_frame.png"

gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.occ.importShapes(r"C:/Users/Liu/Desktop/新建文件夹 (2)/airfoil.stp")
gmsh.model.occ.synchronize()

# 删除不需要的体(recursive: 顺带删掉不再被引用的面/边/点)
gmsh.model.occ.remove([(3, v) for v in REMOVE], recursive=True)
gmsh.model.occ.synchronize()
vols = [t for d, t in gmsh.model.getEntities(3)]
print("remaining volumes:", vols)

# 导出干净 STEP
gmsh.write(OUT_STEP)
print("wrote", OUT_STEP)

# 渲染确认
surf2vol = {}
for v in vols:
    for d, s in gmsh.model.getBoundary([(3, v)], oriented=False):
        surf2vol.setdefault(abs(s), v)
gmsh.option.setNumber("Mesh.MeshSizeMin", 4)
gmsh.option.setNumber("Mesh.MeshSizeMax", 8)
gmsh.model.mesh.generate(2)
ntags, nc, _ = gmsh.model.mesh.getNodes(); nc = nc.reshape(-1, 3)
id2idx = {int(t): i for i, t in enumerate(ntags)}
colors = {1:"#e41a1c",2:"#377eb8",3:"#4daf4a",4:"#984ea3"}
tris = {v: [] for v in vols}
for d, s in gmsh.model.getEntities(2):
    v = surf2vol.get(s)
    if v not in tris: continue
    es = gmsh.model.mesh.getElements(dim=2, tag=s)
    for et, conn in zip(es[0], es[2]):
        npe = gmsh.model.mesh.getElementProperties(et)[3]
        for row in conn.reshape(-1, npe):
            tris[v].append([id2idx[int(n)] for n in row[:3]])

fig = plt.figure(figsize=(16, 8))
for k, (e, a) in enumerate([(20, -60), (20, 120)]):
    ax = fig.add_subplot(1, 2, k+1, projection="3d")
    for v in vols:
        if not tris[v]: continue
        ax.add_collection3d(Poly3DCollection(nc[np.array(tris[v])],
                            facecolor=colors[v], edgecolor=(0,0,0,0.2), linewidths=0.1))
    mn, mx = nc.min(0), nc.max(0)
    ax.set_xlim(mn[0],mx[0]); ax.set_ylim(mn[1],mx[1]); ax.set_zlim(mn[2],mx[2])
    ax.set_box_aspect(mx-mn); ax.view_init(e, a); ax.set_title(f"view {k+1}")
fig.legend(handles=[Patch(color=colors[v], label=f"vol {v}") for v in vols],
           loc="lower center", ncol=4)
plt.tight_layout(rect=[0,0.05,1,1]); plt.savefig(OUT_PNG, dpi=110)
print("tris:", {v: len(tris[v]) for v in vols}, "wrote", OUT_PNG)
gmsh.finalize()
