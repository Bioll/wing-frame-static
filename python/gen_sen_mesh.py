"""
SEN (single-edge notched) tension 试件网格生成器 —— 相场断裂经典 benchmark。

第1步: 均匀结构四边形网格 (CPE4, 平面应变)，含预制缺口 seam(裂纹面节点劈裂)。
导出 Abaqus .inp，自研代码可直接解析。

几何 (单位 mm):
  1x1 方板, [0,1]x[0,1]
  水平预制缺口: 从 (0, 0.5) 到 (0.5, 0.5)，裂纹面上下不共节点
  裂尖在 x=0.5 处闭合(共节点)

节点集:
  TOP    : y=1 上边界
  BOTTOM : y=0 下边界
"""

# ---------------- 参数 ----------------
L      = 1.0      # 板边长 mm
N      = 100      # 每边单元数 (必须偶数, 保证裂尖落在节点上)
CRACK_Y = 0.5     # 缺口高度
CRACK_X = 0.5     # 缺口长度(从左边到此 x)
ELTYPE = "CPE4"   # 平面应变四节点; 平面应力改 CPS4
OUT    = "sen_tension.inp"
# --------------------------------------

assert N % 2 == 0, "N 必须为偶数"
h = L / N
jc = round(CRACK_Y / h)          # 裂纹所在节点行
ic = round(CRACK_X / h)          # 裂尖所在节点列 (x=0.5)

# 节点坐标网格: id = j*(N+1) + i + 1
def base_id(i, j):
    return j * (N + 1) + i + 1

nodes = {}   # id -> (x, y)
for j in range(N + 1):
    for i in range(N + 1):
        nodes[base_id(i, j)] = (i * h, j * h)

# 裂纹面劈裂: 对 j=jc, i=0..ic-1 的节点造一个"上"副本
# 原节点 = 下裂纹面, 副本 = 上裂纹面; i=ic(裂尖)及右侧保持共节点
next_id = (N + 1) * (N + 1) + 1
crack_top = {}   # i -> 上副本节点 id
for i in range(ic):              # i = 0 .. ic-1
    x, y = nodes[base_id(i, jc)]
    crack_top[i] = next_id
    nodes[next_id] = (x, y)
    next_id += 1

# 单元: 行 r 跨节点行 r,r+1; 列 i 跨节点列 i,i+1
# 缺口上方那一行(r=jc)的底边落在裂纹线 -> 开口段用上副本
def node_for(i, j, elem_row):
    # 只有"缺口上方单元(elem_row==jc)的底边节点(j==jc)"且在开口段才换上副本
    if elem_row == jc and j == jc and i in crack_top:
        return crack_top[i]
    return base_id(i, j)

elements = {}    # eid -> (n1,n2,n3,n4)
eid = 1
for r in range(N):
    for i in range(N):
        n1 = node_for(i,   r,   r)
        n2 = node_for(i+1, r,   r)
        n3 = node_for(i+1, r+1, r)
        n4 = node_for(i,   r+1, r)
        elements[eid] = (n1, n2, n3, n4)
        eid += 1

# 边界节点集
top    = [base_id(i, N) for i in range(N + 1)]
bottom = [base_id(i, 0) for i in range(N + 1)]

# ---------------- 写 .inp ----------------
def write_set(f, name, ids):
    f.write(f"*NSET, NSET={name}\n")
    for k in range(0, len(ids), 8):
        f.write(", ".join(str(v) for v in ids[k:k+8]) + ",\n")

with open(OUT, "w") as f:
    f.write("*HEADING\n")
    f.write(f"SEN tension benchmark, uniform {N}x{N} quad, h={h:.4g} mm\n")
    f.write("*NODE\n")
    for nid in sorted(nodes):
        x, y = nodes[nid]
        f.write(f"{nid}, {x:.8f}, {y:.8f}, 0.0\n")
    f.write(f"*ELEMENT, TYPE={ELTYPE}, ELSET=PLATE\n")
    for e in sorted(elements):
        n1, n2, n3, n4 = elements[e]
        f.write(f"{e}, {n1}, {n2}, {n3}, {n4}\n")
    write_set(f, "TOP", top)
    write_set(f, "BOTTOM", bottom)

print(f"h            = {h:.5f} mm")
print(f"裂纹行 jc    = {jc}  (y={jc*h})")
print(f"裂尖列 ic    = {ic}  (x={ic*h})")
print(f"节点数       = {len(nodes)}  (其中裂纹面副本 {len(crack_top)})")
print(f"单元数       = {len(elements)}")
print(f"已写出       = {OUT}")
