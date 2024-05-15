# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import importlib

sys.path.append("/SNS/SNAP/shared/Malcolm/code/crystalBox")
import crystalBox as crys
importlib.reload(crys)

#1. Create a "crystalBox" using:

# a nickname
diamond = crys.Box('diamond')
diamond.tickWS(0.15)

sf = 0.9
diamond.scaleLattice(sf)
print(f"diamond lattice param: {diamond.a:.4f} Ang, vol = {diamond.volume:.4f} Ang^3")


for d in diamond.dSpacing:
    print(d)
diamond.resetLattice()
print(f"diamond lattice param reset to: {diamond.a:.4f} Ang, vol = {diamond.volume:.4f} Ang^3")
for d in diamond.dSpacing:
    print(d)

diamond.plot()