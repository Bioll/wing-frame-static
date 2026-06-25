"""生成悬臂梁四面体网格(结构网格+每六面体劈6四面体), 输出 fem3d 输入格式。
x=0 端固定; 自由端 x=L 加载。LOAD='axial'(轴向,精确解) 或 'bending'(横向)。
"""
import sys

L,B,H = 100.0,10.0,10.0          # 长 x 截面 mm
nx,ny,nz = 20,2,2                # 单元格数
E,NU,RHO = 70000.0,0.3,2.7e-9
P = 1000.0                       # 总载荷 N
LOAD = sys.argv[1] if len(sys.argv)>1 else 'axial'
OUT  = sys.argv[2] if len(sys.argv)>2 else 'cantilever.txt'

def nid(i,j,k): return i + (nx+1)*(j + (ny+1)*k)
nodes=[]
for k in range(nz+1):
    for j in range(ny+1):
        for i in range(nx+1):
            nodes.append((i*L/nx, j*B/ny, k*H/nz))
# 六面体->6四面体(共主对角 0-6)
HEX=[(0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)]
def corners(i,j,k):
    return [nid(i,j,k),nid(i+1,j,k),nid(i+1,j+1,k),nid(i,j+1,k),
            nid(i,j,k+1),nid(i+1,j,k+1),nid(i+1,j+1,k+1),nid(i,j+1,k+1)]
elems=[]
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            c=corners(i,j,k)
            for t in HEX: elems.append((c[t[0]],c[t[1]],c[t[2]],c[t[3]]))

# 固定 x=0 面节点
fix=[n for n,(x,y,z) in enumerate(nodes) if abs(x)<1e-9]
# 自由端 x=L 面节点, 均分载荷
tip=[n for n,(x,y,z) in enumerate(nodes) if abs(x-L)<1e-9]
f_each = P/len(tip)
if LOAD=='axial': fvec=(f_each,0,0)
else:             fvec=(0,-f_each,0)   # bending: -y

with open(OUT,'w') as f:
    f.write("NODES %d\n"%len(nodes))
    for x,y,z in nodes: f.write("%.10g %.10g %.10g\n"%(x,y,z))
    f.write("ELEMENTS %d\n"%len(elems))
    for e in elems: f.write("%d %d %d %d\n"%e)
    f.write("MATERIAL %g %g %g\n"%(E,NU,RHO))
    f.write("GRAVITY 0 0 0\n")
    f.write("FIX %d\n"%len(fix))
    for n in fix: f.write("%d 1 1 1\n"%n)
    f.write("NFORCE %d\n"%len(tip))
    for n in tip: f.write("%d %g %g %g\n"%(n,fvec[0],fvec[1],fvec[2]))

# 理论解
A=B*H; I=B*H**3/12
print("LOAD=%s nodes=%d elems=%d"%(LOAD,len(nodes),len(elems)))
if LOAD=='axial':
    print("理论轴向伸长 PL/(AE) = %.6g mm"%(P*L/(A*E)))
else:
    print("理论梁端挠度 PL^3/(3EI) = %.6g mm"%(P*L**3/(3*E*I)))
