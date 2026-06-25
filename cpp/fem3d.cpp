// fem3d.cpp —— 最小线弹性有限元求解器 (4节点四面体 C3D4)
// 组装 -> Dirichlet 约束 -> 共轭梯度求解 -> 输出 VTK(位移+von Mises)
// 自包含, 无外部依赖。编译: g++ -O2 -std=c++17 fem3d.cpp -o fem3d
//
// 输入格式(文本):
//   NODES n / 每行: x y z
//   ELEMENTS m / 每行: n0 n1 n2 n3   (0基节点号)
//   MATERIAL E nu rho
//   GRAVITY gx gy gz                 (重力加速度矢量, 体力)
//   FIX k / 每行: node cx cy cz       (c=1 表示该自由度固定)
//   NFORCE p / 每行: node fx fy fz    (集中力, 可选)
// 用法: fem3d input.txt output.vtk
#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x,y,z; };

// 3x3 行列式与逆
static double det3(const double m[3][3]){
    return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
         - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
         + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
}
static void inv3(const double m[3][3], double inv[3][3], double det){
    double id=1.0/det;
    inv[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])*id;
    inv[0][1]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])*id;
    inv[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])*id;
    inv[1][0]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])*id;
    inv[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])*id;
    inv[1][2]=(m[0][2]*m[1][0]-m[0][0]*m[1][2])*id;
    inv[2][0]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])*id;
    inv[2][1]=(m[0][1]*m[2][0]-m[0][0]*m[2][1])*id;
    inv[2][2]=(m[0][0]*m[1][1]-m[0][1]*m[1][0])*id;
}

int main(int argc,char**argv){
    if(argc<3){ fprintf(stderr,"用法: fem3d input.txt output.vtk\n"); return 1; }
    ifstream in(argv[1]);
    string tok;
    int nn=0, ne=0;
    vector<Vec3> X; vector<array<int,4>> E;
    double Em=0,nu=0,rho=0; double g[3]={0,0,0};
    vector<array<int,4>> fixspec; vector<array<double,4>> nforce;
    while(in>>tok){
        if(tok=="NODES"){ in>>nn; X.resize(nn);
            for(int i=0;i<nn;i++) in>>X[i].x>>X[i].y>>X[i].z; }
        else if(tok=="ELEMENTS"){ in>>ne; E.resize(ne);
            for(int i=0;i<ne;i++) in>>E[i][0]>>E[i][1]>>E[i][2]>>E[i][3]; }
        else if(tok=="MATERIAL"){ in>>Em>>nu>>rho; }
        else if(tok=="GRAVITY"){ in>>g[0]>>g[1]>>g[2]; }
        else if(tok=="FIX"){ int k; in>>k; for(int i=0;i<k;i++){ array<int,4> a; in>>a[0]>>a[1]>>a[2]>>a[3]; fixspec.push_back(a);} }
        else if(tok=="NFORCE"){ int p; in>>p; for(int i=0;i<p;i++){ array<double,4> a; in>>a[0]>>a[1]>>a[2]>>a[3]; nforce.push_back(a);} }
    }
    int ndof=3*nn;
    printf("nodes=%d elems=%d dof=%d\n",nn,ne,ndof);

    // 本构 D (6x6), Voigt: [xx,yy,zz,xy,yz,zx]
    double D[6][6]={0};
    double c=Em/((1+nu)*(1-2*nu));
    D[0][0]=D[1][1]=D[2][2]=c*(1-nu);
    D[0][1]=D[0][2]=D[1][0]=D[1][2]=D[2][0]=D[2][1]=c*nu;
    D[3][3]=D[4][4]=D[5][5]=c*(1-2*nu)/2;

    // 稀疏K(三元组->map)与载荷
    map<long long,double> K;
    vector<double> F(ndof,0.0);
    auto addK=[&](int i,int j,double v){ K[(long long)i*ndof+j]+=v; };

    for(int e=0;e<ne;e++){
        int id[4]={E[e][0],E[e][1],E[e][2],E[e][3]};
        Vec3 p0=X[id[0]],p1=X[id[1]],p2=X[id[2]],p3=X[id[3]];
        double J[3][3]={ {p1.x-p0.x,p2.x-p0.x,p3.x-p0.x},
                         {p1.y-p0.y,p2.y-p0.y,p3.y-p0.y},
                         {p1.z-p0.z,p2.z-p0.z,p3.z-p0.z} };
        double dJ=det3(J); double V=dJ/6.0;
        if(V<0){ /* 节点序反向, 体积取正 */ }
        double Ji[3][3]; inv3(J,Ji,dJ);
        // dN/dξ (4x3)
        double dNr[4][3]={ {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} };
        // dN/dx = dNr * Ji
        double dN[4][3];
        for(int a=0;a<4;a++)for(int j=0;j<3;j++){ double s=0; for(int b=0;b<3;b++) s+=dNr[a][b]*Ji[b][j]; dN[a][j]=s; }
        // B (6x12)
        double B[6][12]={0};
        for(int a=0;a<4;a++){
            double bx=dN[a][0],by=dN[a][1],bz=dN[a][2];
            int c0=3*a;
            B[0][c0+0]=bx; B[1][c0+1]=by; B[2][c0+2]=bz;
            B[3][c0+0]=by; B[3][c0+1]=bx;
            B[4][c0+1]=bz; B[4][c0+2]=by;
            B[5][c0+0]=bz; B[5][c0+2]=bx;
        }
        double vol=fabs(V);
        // Ke = vol * B^T D B
        double DB[6][12];
        for(int r=0;r<6;r++)for(int cc=0;cc<12;cc++){ double s=0; for(int k=0;k<6;k++) s+=D[r][k]*B[k][cc]; DB[r][cc]=s; }
        for(int a=0;a<12;a++)for(int b=0;b<12;b++){ double s=0; for(int k=0;k<6;k++) s+=B[k][a]*DB[k][b]; double v=s*vol;
            int gi=3*id[a/3]+a%3, gj=3*id[b/3]+b%3; addK(gi,gj,v); }
        // 体力(重力): 一致节点力 rho*g*vol/4
        for(int a=0;a<4;a++) for(int d=0;d<3;d++) F[3*id[a]+d]+=rho*g[d]*vol/4.0;
    }
    // 集中力
    for(auto&a:nforce){ int n=(int)a[0]; F[3*n+0]+=a[1]; F[3*n+1]+=a[2]; F[3*n+2]+=a[3]; }

    // Dirichlet: 固定自由度置零行列, 对角1, F=0 (保持对称)
    vector<char> fixed(ndof,0);
    for(auto&a:fixspec){ int n=a[0]; for(int d=0;d<3;d++) if(a[1+d]) fixed[3*n+d]=1; }
    // 先把固定列对自由侧的贡献移到右端(此处规定位移=0, 无需移项), 然后清行列
    for(auto it=K.begin();it!=K.end();){
        int i=it->first/ndof, j=it->first%ndof;
        if(fixed[i]||fixed[j]){ it=K.erase(it); } else ++it;
    }
    for(int i=0;i<ndof;i++) if(fixed[i]){ K[(long long)i*ndof+i]=1.0; F[i]=0.0; }

    // 转 CSR
    vector<int> rp(ndof+1,0); vector<int> ci; vector<double> va;
    {
        vector<int> cnt(ndof,0);
        for(auto&kv:K) cnt[kv.first/ndof]++;
        for(int i=0;i<ndof;i++) rp[i+1]=rp[i]+cnt[i];
        ci.resize(rp[ndof]); va.resize(rp[ndof]);
        vector<int> pos(rp.begin(),rp.end()-1);
        for(auto&kv:K){ int i=kv.first/ndof,j=kv.first%ndof; ci[pos[i]]=j; va[pos[i]]=kv.second; pos[i]++; }
    }
    auto matvec=[&](const vector<double>&x, vector<double>&y){
        for(int i=0;i<ndof;i++){ double s=0; for(int k=rp[i];k<rp[i+1];k++) s+=va[k]*x[ci[k]]; y[i]=s; }
    };
    // Jacobi 预处理对角
    vector<double> diag(ndof,1.0);
    for(int i=0;i<ndof;i++) for(int k=rp[i];k<rp[i+1];k++) if(ci[k]==i) diag[i]=va[k];

    // 预处理共轭梯度
    vector<double> u(ndof,0), r=F, z(ndof), p(ndof), Ap(ndof);
    for(int i=0;i<ndof;i++) z[i]=r[i]/diag[i];
    p=z;
    double rz=0; for(int i=0;i<ndof;i++) rz+=r[i]*z[i];
    double bnorm=0; for(double v:F) bnorm+=v*v; bnorm=sqrt(bnorm)+1e-30;
    int it=0, maxit=20000; double tol=1e-8;
    for(;it<maxit;it++){
        matvec(p,Ap);
        double pAp=0; for(int i=0;i<ndof;i++) pAp+=p[i]*Ap[i];
        double al=rz/pAp;
        double rn=0;
        for(int i=0;i<ndof;i++){ u[i]+=al*p[i]; r[i]-=al*Ap[i]; rn+=r[i]*r[i]; }
        if(sqrt(rn)/bnorm<tol) break;
        for(int i=0;i<ndof;i++) z[i]=r[i]/diag[i];
        double rz2=0; for(int i=0;i<ndof;i++) rz2+=r[i]*z[i];
        double be=rz2/rz; rz=rz2;
        for(int i=0;i<ndof;i++) p[i]=z[i]+be*p[i];
    }
    printf("CG 迭代=%d, 相对残差=%.2e\n",it,sqrt([&]{double s=0;for(double v:r)s+=v*v;return s;}())/bnorm);

    // 单元 von Mises
    vector<double> vm(ne,0.0), umag(nn,0.0);
    for(int i=0;i<nn;i++) umag[i]=sqrt(u[3*i]*u[3*i]+u[3*i+1]*u[3*i+1]+u[3*i+2]*u[3*i+2]);
    double maxu=0,maxvm=0;
    for(int e=0;e<ne;e++){
        int id[4]={E[e][0],E[e][1],E[e][2],E[e][3]};
        Vec3 p0=X[id[0]],p1=X[id[1]],p2=X[id[2]],p3=X[id[3]];
        double J[3][3]={ {p1.x-p0.x,p2.x-p0.x,p3.x-p0.x},{p1.y-p0.y,p2.y-p0.y,p3.y-p0.y},{p1.z-p0.z,p2.z-p0.z,p3.z-p0.z} };
        double dJ=det3(J); double Ji[3][3]; inv3(J,Ji,dJ);
        double dNr[4][3]={ {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} }; double dN[4][3];
        for(int a=0;a<4;a++)for(int j=0;j<3;j++){ double s=0; for(int b=0;b<3;b++) s+=dNr[a][b]*Ji[b][j]; dN[a][j]=s; }
        double B[6][12]={0};
        for(int a=0;a<4;a++){ double bx=dN[a][0],by=dN[a][1],bz=dN[a][2]; int c0=3*a;
            B[0][c0]=bx;B[1][c0+1]=by;B[2][c0+2]=bz;B[3][c0]=by;B[3][c0+1]=bx;B[4][c0+1]=bz;B[4][c0+2]=by;B[5][c0]=bz;B[5][c0+2]=bx; }
        double ue[12]; for(int a=0;a<4;a++)for(int d=0;d<3;d++) ue[3*a+d]=u[3*id[a]+d];
        double eps[6]={0}; for(int r=0;r<6;r++)for(int cc=0;cc<12;cc++) eps[r]+=B[r][cc]*ue[cc];
        double s[6]={0}; for(int r=0;r<6;r++)for(int k=0;k<6;k++) s[r]+=D[r][k]*eps[k];
        double mis=sqrt(0.5*((s[0]-s[1])*(s[0]-s[1])+(s[1]-s[2])*(s[1]-s[2])+(s[2]-s[0])*(s[2]-s[0]))+3*(s[3]*s[3]+s[4]*s[4]+s[5]*s[5]));
        vm[e]=mis; maxvm=max(maxvm,mis);
    }
    for(int i=0;i<nn;i++) maxu=max(maxu,umag[i]);
    printf("最大|U|=%.6g  最大vonMises=%.6g\n",maxu,maxvm);

    // 输出 VTK
    FILE*f=fopen(argv[2],"w");
    fprintf(f,"# vtk DataFile Version 3.0\nfem3d\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(f,"POINTS %d float\n",nn);
    for(int i=0;i<nn;i++) fprintf(f,"%g %g %g\n",X[i].x,X[i].y,X[i].z);
    fprintf(f,"CELLS %d %d\n",ne,ne*5);
    for(int e=0;e<ne;e++) fprintf(f,"4 %d %d %d %d\n",E[e][0],E[e][1],E[e][2],E[e][3]);
    fprintf(f,"CELL_TYPES %d\n",ne); for(int e=0;e<ne;e++) fprintf(f,"10\n");
    fprintf(f,"POINT_DATA %d\nVECTORS U float\n",nn);
    for(int i=0;i<nn;i++) fprintf(f,"%g %g %g\n",u[3*i],u[3*i+1],u[3*i+2]);
    fprintf(f,"SCALARS Umag float 1\nLOOKUP_TABLE default\n");
    for(int i=0;i<nn;i++) fprintf(f,"%g\n",umag[i]);
    fprintf(f,"CELL_DATA %d\nSCALARS Mises float 1\nLOOKUP_TABLE default\n",ne);
    for(int e=0;e<ne;e++) fprintf(f,"%g\n",vm[e]);
    fclose(f);
    printf("已写出 %s\n",argv[2]);
    return 0;
}
