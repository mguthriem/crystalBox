# crystalFOM: a utility to calculate a figure of merit, reflecting "measurability" of a powder sample on SNAP
# M. Guthrie 20241204

import numpy as np
import sys
sys.path.append("/SNS/SNAP/shared/Malcolm/code/crystalBox/")
import crystalBox as crys
import importlib
importlib.reload(crys)

cifList = ['diamond','nickel', 'kdp', 'dkdp']

print("CIF source, FOM, abs/inc FOM")
for cif in cifList:
    
    crystal = crys.Box(cif)
    fom1,fom2 = crystal.calcFOM()
    print(f"{cif}, {fom1:.4f} {fom2:.4f}")