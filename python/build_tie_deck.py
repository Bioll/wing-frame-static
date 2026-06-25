"""由非合并四面体网格构建 TIE 绑定的 Abaqus 静力 deck。
- 连通块分件 -> 4 个 ELSET
- 件间重合界面 -> element-based surface -> *TIE (面-面绑定, 带位置容差)
- 铝合金 / 翼根固定 / 重力+顶面压力
单位制: mm-N-MPa-tonne-s
"""
import numpy as np
from collections import defaultdict, Counter

SRC="airfoil_parts_tet.inp"; DST="airfoil_tie.inp"
# EDIT 物理参数
E,NU,RHO=70000.0,0.33,2.7e-9
G_ACC=9810.0; G_DIR=(0.0,-1.0,0.0)
PRESSURE=0.01; ROOT_Z=40.0; TOP_NY=0.3
TOL=1e-3

# ---- 读网格 ----
nodes={}; elems=[]; sec=None
for line in open(SRC,encoding="utf-8",errors="ignore"):
    t=line.strip()
    if t.startswith("*"):
        u=t.upper(); sec="N" if u.startswith("*NODE") else ("E" if u.startswith("*ELEMENT") else None); continue
    if not t: continue
    p=t.replace(","," ").split()
    if sec=="N":
        try: nodes[int(p[0])]=np.array([float(p[1]),float(p[2]),float(p[3])])
        except: pass
    elif sec=="E":
        try:
            q=[int(x) for x in p]
            if len(q)==5: elems.append(q)
        except: pass

# ---- 连通块分件 ----
parent={}
def find(x):
    parent.setdefault(x,x)
    while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
    return x
def uni(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: parent[ra]=rb
for e in elems:
    for n in e[2:]: uni(e[1],n)
root_of={};
labels={}
order=[r for r,_ in Counter(find(e[1]) for e in elems).most_common()]
lab={r:i+1 for i,r in enumerate(order)}
part_of_elem={e[0]:lab[find(e[1])] for e in elems}
parts=sorted(set(part_of_elem.values()))
print("件数:",len(parts),{p:sum(1 for v in part_of_elem.values() if v==p) for p in parts})

# ---- 每件边界面 (件内只出现1次的面) ----
# C3D4 面: S1=123 S2=124 S3=234 S4=134; 存(对面节点用于法向定向)
FACES=[(1,2,3),(1,2,4),(2,3,4),(1,3,4)]; OPP=[4,3,1,2]
econn={e[0]:e[1:] for e in elems}
face_in_part=defaultdict(list)   # (part, nodekey) -> [(eid,sidx)]
for e in elems:
    eid=e[0]; nd=e[1:]; pt=part_of_elem[eid]
    for si,(a,b,c) in enumerate(FACES,1):
        key=tuple(sorted((nd[a-1],nd[b-1],nd[c-1])))
        face_in_part[(pt,key)].append((eid,si))
bnd=defaultdict(list)            # part -> [(eid,si,nodekey)]
for (pt,key),lst in face_in_part.items():
    if len(lst)==1:
        eid,si=lst[0]; bnd[pt].append((eid,si,key))

# ---- 件间界面: 按邻近(非共形) ----
TOL_IF=1.8
def fcentroid(eid,si):
    a,b,c=FACES[si-1]; nd=econn[eid]
    return (nodes[nd[a-1]]+nodes[nd[b-1]]+nodes[nd[c-1]])/3.0
bnode_xyz={}
for pt,lst in bnd.items():
    ns=set()
    for eid,si,key in lst: ns.update(key)
    bnode_xyz[pt]=np.array([nodes[n] for n in ns])
def cents_of(pt):
    lst=bnd[pt]
    return np.array([fcentroid(e,s) for e,s,k in lst]), [(e,s) for e,s,k in lst]
def near_mask(cents, ref):
    out=np.zeros(len(cents),bool)
    for k in range(0,len(cents),200):
        ch=cents[k:k+200]
        d=np.sqrt(((ch[:,None,:]-ref[None,:,:])**2).sum(-1)).min(1)
        out[k:k+200]=d<TOL_IF
    return out
iface={}; interface_faces=set()
for i in parts:
    ci,fi=cents_of(i)
    for j in parts:
        if j<=i: continue
        cj,fj=cents_of(j)
        mi=near_mask(ci,bnode_xyz[j]); mj=near_mask(cj,bnode_xyz[i])
        fi_if=[fi[k] for k in range(len(fi)) if mi[k]]
        fj_if=[fj[k] for k in range(len(fj)) if mj[k]]
        if fi_if and fj_if:
            iface[(i,j)]=(fi_if,fj_if)
            interface_faces.update(fi_if); interface_faces.update(fj_if)
pairs=sorted(iface)
print("TIE界面对:",[(k,len(iface[k][0]),len(iface[k][1])) for k in pairs])

# ---- 顶面(外表面非界面, 法向+Y) 加压 ----
top=[]
for pt,lst in bnd.items():
    for eid,si,key in lst:
        if (eid,si) in interface_faces: continue
        a,b,c=FACES[si-1]; nd=econn[eid]
        v1,v2,v3=nodes[nd[a-1]],nodes[nd[b-1]],nodes[nd[c-1]]; vo=nodes[nd[OPP[si-1]-1]]
        n=np.cross(v2-v1,v3-v1)
        if np.dot(n,v1-vo)<0: n=-n
        if n[1]/(np.linalg.norm(n)+1e-12) > TOP_NY: top.append((eid,si))

# 翼根节点
root=[nid for nid,xyz in nodes.items() if xyz[2]<ROOT_Z]
print("root nodes:",len(root)," top faces:",len(top))

# ---- 写 deck ----
def wsurf(f,name,faces):
    bysn=defaultdict(list)
    for eid,si in faces: bysn[si].append(eid)
    for si,els in sorted(bysn.items()):
        f.write(f"*ELSET, ELSET={name}_S{si}\n")
        for i in range(0,len(els),10): f.write(", ".join(map(str,els[i:i+10]))+"\n")
    f.write(f"*SURFACE, TYPE=ELEMENT, NAME={name}\n")
    for si in sorted(bysn): f.write(f"{name}_S{si}, S{si}\n")

with open(DST,"w") as f:
    f.write("*HEADING\n airfoil frame static, parts TIE-bonded, aluminum\n")
    f.write("*NODE\n")
    for nid in sorted(nodes):
        x,y,z=nodes[nid]; f.write(f"{nid}, {x:.8f}, {y:.8f}, {z:.8f}\n")
    # 各件 ELEMENT/ELSET
    by_part=defaultdict(list)
    for e in elems: by_part[part_of_elem[e[0]]].append(e)
    for pt in parts:
        f.write(f"*ELEMENT, TYPE=C3D4, ELSET=PART{pt}\n")
        for e in by_part[pt]: f.write(", ".join(map(str,e))+"\n")
    f.write("*ELSET, ELSET=ALL_FRAME\n")
    f.write(", ".join(f"PART{pt}" for pt in parts)+"\n")
    # 界面 surface + TIE: 只连到中枢件(parts[0]), 中枢当主面, 其他当从面
    # -> 每个从节点只属于一个 tie, 避免过约束
    hub=parts[0]
    for (pa,pb) in pairs:
        if hub not in (pa,pb):
            print(f"  跳过非中枢TIE {pa}-{pb}(已通过件{hub}连通)"); continue
        other = pb if pa==hub else pa
        fhub = iface[(pa,pb)][0] if pa==hub else iface[(pa,pb)][1]
        foth = iface[(pa,pb)][1] if pa==hub else iface[(pa,pb)][0]
        mn,sn=f"IF_{hub}_{other}_MST",f"IF_{hub}_{other}_SLV"
        wsurf(f,mn,sorted(set(fhub))); wsurf(f,sn,sorted(set(foth)))
        f.write(f"*TIE, NAME=TIE_{hub}_{other}, POSITION TOLERANCE=1.0, ADJUST=NO\n{sn}, {mn}\n")
    # 材料/截面
    f.write("*MATERIAL, NAME=ALUMINUM\n*ELASTIC\n")
    f.write(f"{E}, {NU}\n*DENSITY\n{RHO}\n")
    f.write("*SOLID SECTION, ELSET=ALL_FRAME, MATERIAL=ALUMINUM\n")
    # 翼根固定
    f.write("*NSET, NSET=ROOT\n")
    for i in range(0,len(root),10): f.write(", ".join(map(str,root[i:i+10]))+"\n")
    f.write("*BOUNDARY\nROOT, ENCASTRE\n")
    # 顶面压力面
    wsurf(f,"TOP_SURF",top)
    # 静力步 (无 STABILIZE, 靠 TIE 连接)
    f.write("*STEP\n*STATIC\n")
    gx,gy,gz=G_DIR
    f.write(f"*DLOAD\nALL_FRAME, GRAV, {G_ACC}, {gx}, {gy}, {gz}\n")
    f.write(f"*DSLOAD\nTOP_SURF, P, {PRESSURE}\n")
    f.write("*OUTPUT, FIELD\n*NODE OUTPUT\nU, RF\n*ELEMENT OUTPUT\nS, E\n*END STEP\n")
print("wrote",DST)
