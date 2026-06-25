from odbAccess import openOdb
import numpy as np
o=openOdb('airfoil_tie.odb')
fr=o.steps['Step-1'].frames[-1]
U=fr.fieldOutputs['U']; S=fr.fieldOutputs['S']; RF=fr.fieldOutputs['RF']
um=max(v.magnitude for v in U.values)
mis=max(v.mises for v in S.values)
rf=np.array([v.data for v in RF.values]); tot=rf.sum(0)
print('MAXU_mm=%.4g' % um)
print('MAXMISES_MPa=%.4g' % mis)
print('ROOT_RF_N=%.3g,%.3g,%.3g' % (tot[0],tot[1],tot[2]))
o.close()
