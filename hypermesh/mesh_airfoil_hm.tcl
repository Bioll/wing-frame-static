# =====================================================================
# HyperMesh 2021: 读 airfoil.stp -> 面网格 -> 导出 Abaqus .inp
# 命令签名按 install 自带 aerospace/SurfaceMesh 工具校准, hmbatch 实测可跑。
# GUI 运行: File -> Run -> Tcl/Tk Script
# 批处理:   hmbatch.exe -tcl mesh_airfoil_hm.tcl
# =====================================================================
set f    "C:/Users/Liu/Desktop/新建文件夹 (2)/airfoil.stp"
set out  "C:/Users/Liu/desktop/SEN_phasefield/airfoil_hm.inp"
set tpl  "C:/Program Files/Altair/2021/hwdesktop/templates/feoutput/abaqus/standard.3d"
set size 8.0

# 1. 导入几何 (STEP/STEP_CT 翻译器在本机坏了, 用 AutoDetect)
*geomimport "AutoDetect" "$f"
*createmark surfaces 1 "all"
set surfList [hm_getmark surfaces 1]
set numSurfs [llength $surfList]
puts "imported surfaces = $numSurfs"

# 2. 对所有曲面 automesh (size mm, 三角单元)
*setgeomrefinelevel 1
*setedgedensitylinkwithaspectratio -1
*elementorder 1
eval *createmark surfaces 1 $surfList
*interactiveremeshsurf 1 $size 2 2 2 1 1
for {set i 0} {$i < $numSurfs} {incr i} {
    *set_meshfaceparams $i 2 2 0 0 1 0.5 1 1
    *automesh $i 2 2
}
*storemeshtodatabase 1
*ameshclearsurface
*createmark elements 1 "all"
puts "elements = [hm_marklength elements 1]"

# 3. 导出 Abaqus .inp
*createmark elements 1 "all"
*feoutputwithdata "$tpl" "$out" 0 0 1 1 0
puts "exported: $out  exists=[file exists $out]"
puts "DONE"
