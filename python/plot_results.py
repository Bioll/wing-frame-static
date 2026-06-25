# abaqus viewer noGUI=plot_results.py  -> 出位移/应力云图 PNG
from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb
import visualization

odb = session.openOdb(name='airfoil_tie.odb')
vp = session.viewports['Viewport: 1']
vp.setValues(displayedObject=odb, width=240, height=150)
vp.viewportAnnotationOptions.setValues(triad=OFF, legendDecimalPlaces=3,
                                       compass=OFF, title=OFF, state=OFF)
vp.odbDisplay.commonOptions.setValues(visibleEdges=FREE)
session.pngOptions.setValues(imageSize=(1700, 1050))
session.printOptions.setValues(vpDecorations=OFF, reduceColors=False)

# 标准等轴视图
vp.view.setValues(session.views['Iso'])

# 1) 位移幅值(变形图)
vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF,))
vp.odbDisplay.setPrimaryVariable(variableLabel='U', outputPosition=NODAL,
                                 refinement=(INVARIANT, 'Magnitude'))
vp.view.fitView()
session.printToFile(fileName='res_disp', format=PNG, canvasObjects=(vp,))

# 2) von Mises 应力
vp.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT,
                                 refinement=(INVARIANT, 'Mises'))
vp.view.fitView()
session.printToFile(fileName='res_mises', format=PNG, canvasObjects=(vp,))

# 最大值
fr = odb.steps['Step-1'].frames[-1]
um = max(v.magnitude for v in fr.fieldOutputs['U'].values)
mis = max(v.mises for v in fr.fieldOutputs['S'].values)
print('RESULT MAXU_mm=%.4g MAXMISES_MPa=%.4g' % (um, mis))
odb.close()
