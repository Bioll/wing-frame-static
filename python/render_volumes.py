"""把 7 个体按颜色分开渲染, 便于辨认蒙皮。argv 给要隐藏的体号。"""
import sys, gmsh, numpy as np
HIDE = set(int(x) for x in sys.argv[1:])
OUTPNG = "airfoil_volumes.png" if not HIDE else f"airfoil_hide{'_'.join(map(str,sorted(HIDE)))}.png"
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Patch

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.occ.importShapes(r"C:/Users/Liu/Desktop/新建文件夹 (2)/airfoil.stp")
gmsh.model.occ.synchronize()

# surface -> volume 映射(共享面取第一个体)
surf2vol = {}
vols = [t for d, t in gmsh.model.getEntities(3)]
for v in vols:
    for d, s in gmsh.model.getBoundary([(3, v)], oriented=False):
        surf2vol.setdefault(abs(s), v)

gmsh.option.setNumber("Mesh.MeshSizeMin", 6)
gmsh.option.setNumber("Mesh.MeshSizeMax", 12)
gmsh.model.mesh.generate(2)

ntags, ncoords, _ = gmsh.model.mesh.getNodes()
ncoords = ncoords.reshape(-1, 3)
id2idx = {int(t): i for i, t in enumerate(ntags)}

colors = {1:"#e41a1c",2:"#377eb8",3:"#4daf4a",4:"#984ea3",
          5:"#ff7f00",6:"#a65628",7:"#f0e000"}
tris_by_vol = {v: [] for v in vols}
for d, s in gmsh.model.getEntities(2):
    v = surf2vol.get(s)
    if v is None: continue
    ets, etg, enod = gmsh.model.mesh.getElements(dim=2, tag=s)
    for et, conn in zip(ets, enod):
        npe = gmsh.model.mesh.getElementProperties(et)[3]
        conn = conn.reshape(-1, npe)
        for row in conn:
            tris_by_vol[v].append([id2idx[int(n)] for n in row[:3]])
gmsh.finalize()

fig = plt.figure(figsize=(16, 8))
for k, (elev, azim) in enumerate([(20, -60), (20, 120)]):
    ax = fig.add_subplot(1, 2, k+1, projection="3d")
    for v in vols:
        if not tris_by_vol[v] or v in HIDE: continue
        verts = ncoords[np.array(tris_by_vol[v])]
        ax.add_collection3d(Poly3DCollection(
            verts, facecolor=colors[v], edgecolor=(0,0,0,0.25),
            linewidths=0.1, alpha=0.95))
    mn, mx = ncoords.min(0), ncoords.max(0)
    ax.set_xlim(mn[0],mx[0]); ax.set_ylim(mn[1],mx[1]); ax.set_zlim(mn[2],mx[2])
    ax.set_box_aspect(mx-mn); ax.view_init(elev=elev, azim=azim)
    ax.set_title(f"view {k+1}")
leg = [Patch(color=colors[v], label=f"vol {v}") for v in vols]
fig.legend(handles=leg, loc="lower center", ncol=7)
plt.tight_layout(rect=[0,0.05,1,1])
plt.savefig(OUTPNG, dpi=110)
print("hidden:", sorted(HIDE), "counts:", {v: len(tris_by_vol[v]) for v in vols})
print("wrote", OUTPNG)
