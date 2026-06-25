# 不合并界面的整体四面体(去掉 fix_comp_bdr) -> 4件天然成4个不相连连通块, Python再分件
set f    "C:/Users/Liu/desktop/SEN_phasefield/airfoil_frame.stp"
set out  "C:/Users/Liu/desktop/SEN_phasefield/airfoil_parts_tet.inp"
set tpl  "C:/Program Files/Altair/2021/hwdesktop/templates/feoutput/abaqus/standard.3d"
set size 8.0
set log  [open "C:/Users/Liu/desktop/SEN_phasefield/sep_log.txt" w]
proc L {fp m} { puts $fp $m; flush $fp }

*geomimport "AutoDetect" "$f"
*createmark surfaces 1 "all"
set surfList [hm_getmark surfaces 1]
set numSurfs [llength $surfList]
*setgeomrefinelevel 1
*setedgedensitylinkwithaspectratio -1
*elementorder 1
eval *createmark surfaces 1 $surfList
*interactiveremeshsurf 1 $size 0 0 2 1 1
for {set i 0} {$i < $numSurfs} {incr i} {
    *set_meshfaceparams $i 0 0 0 0 1 0.5 1 1
    *automesh $i 0 0
}
*storemeshtodatabase 1
*ameshclearsurface
*createmark elements 1 "by config" 103
L $log "shellTris=[hm_marklength elements 1]"

# 结构 tetra 不含 fix_comp_bdr -> 界面不合并
*createstringarray 2 "pars: upd_shell tet_clps='0.300000,0.600000,0.900000,1.000000'" "tet: 35 1.2 -1 0 0.8 0"
*createmark components 1 all
if {[catch {*tetmesh components 1 9 elements 0 -1 1 2} e]} { L $log "tet ERR: $e" }
*createmark elements 1 "by config" 204
L $log "tet3D=[hm_marklength elements 1]"
*createmark elements 1 "by config" 103
if {[hm_marklength elements 1] > 0} { *deletemark elements 1 }
*createmark elements 1 "all"
*feoutputwithdata "$tpl" "$out" 0 0 1 1 0
L $log "exported=[file exists $out]"
close $log
