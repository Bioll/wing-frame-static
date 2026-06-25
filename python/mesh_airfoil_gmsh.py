"""
用 gmsh 给 airfoil.stp 划网格 + 渲染出图 + 导出 .inp。
第1步先做表面网格(三角形)用于可视化, 体网格(四面体)留作下一步。
"""
import gmsh, math
import matplotlib
matplotlib.use("Agg")          # 无界面后端, 直接出 png
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

STEP = r"C:/Users/Liu/Desktop/新建文件夹 (2)/airfoil.stp"
CLEN = 8.0                     # 目标单元尺寸 mm (出图用, 可调)
OUT_INP = "airfoil_surf.inp"
OUT_PNG = "airfoil_mesh.png"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.occ.importShapes(STEP)
gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.MeshSizeMin", CLEN * 0.5)
gmsh.option.setNumber("Mesh.MeshSizeMax", CLEN)
gmsh.model.mesh.generate(2)    # 2D = 表面网格

# 导出 inp (Abaqus)
gmsh.write(OUT_INP)

# 取出三角形用于渲染
ntags, ncoords, _ = gmsh.model.mesh.getNodes()
ncoords = ncoords.reshape(-1, 3)
id2idx = {int(t): i for i, t in enumerate(ntags)}
etypes, etags, enodes = gmsh.model.mesh.getElements(dim=2)
tris = []
for et, conn in zip(etypes, enodes):
    npe = gmsh.model.mesh.getElementProperties(et)[3]
    conn = conn.reshape(-1, npe)
    for row in conn:
        tris.append([id2idx[int(n)] for n in row[:3]])
tris = np.array(tris)
print("nodes:", len(ntags), " surface elems:", len(tris))
gmsh.finalize()

# 渲染
fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111, projection="3d")
verts = ncoords[tris]
pc = Poly3DCollection(verts, facecolor=(0.6, 0.75, 0.9),
                      edgecolor=(0.15, 0.15, 0.15), linewidths=0.15)
ax.add_collection3d(pc)
xmin, ymin, zmin = ncoords.min(0); xmax, ymax, zmax = ncoords.max(0)
ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax); ax.set_zlim(zmin, zmax)
ax.set_box_aspect((xmax-xmin, ymax-ymin, zmax-zmin))
ax.view_init(elev=20, azim=-60)
ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
ax.set_title(f"airfoil surface mesh (gmsh, clen={CLEN}mm, {len(tris)} tris)")
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=110)
print("wrote", OUT_PNG, "and", OUT_INP)
