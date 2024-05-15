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
nickel = crys.Box('nickel')

# or if the cif file lives in /SNS/SNAP/shared/cifLibrary, just the file name works
NAC = crys.Box('Na2Ca3Al2F14.cif')

# or a full path to a cif file
NAC = crys.Box('/SNS/SNAP/shared/cifLibrary/Na2Ca3Al2F14.cif')

#2. Display summary of crystal information using

print("A summary of the diamond crystal info:")
print(diamond.summary())

#3. get a list of useful numbers, such as d-spacing

print("\nThese are d-spacings for nickel:")
for i in range(nickel.nRef):
    print(f"{nickel.dSpacing[i]:.4f}")

#4. you also have these other useful numbers:

print("\nYou can also get: hkl, d-spacing multiplicity and estimated intensity. e.g. for NAC:")
for i in range(10):
    print(f"{i} {NAC.hkl[i]} {NAC.dSpacing[i]:.4f} {NAC.fSq[i]:.1f} {NAC.mult[i]} {NAC.estInt[i]:.1f}")
    
#5. You can apply a simple linear scaling to the lattice parameters using

print(f"\noriginal nickel lattice param: {nickel.a:.4f} Ang, vol = {nickel.volume:.4f} Ang^3")
sf = 0.75
nickel.scaleLattice(sf)
print(f"nickel lattice param scaled by {sf}: {nickel.a:.4f} Ang, vol = {nickel.volume:.4f} Ang^3")

nickel.resetLattice()

print(f"reset nickel lattice param: {nickel.a:.4f} Ang, vol = {nickel.volume:.4f} Ang^3")
print(nickel.dSpacing)

#6. You can create a workspace with ticks using

print(f"\nNow creating workspace \'ticks: diamond\' with diamond peak positions")
diamond.tickWS(-0.01) #the value in brackets is the height of the tick marks

#and then plot is using

diamond.plot()

#you can then drag diffraction data onto that plot

LoadNexus(Filename='/SNS/SNAP/shared/Malcolm/code/crystalBox/data/diamond.nxs',
    OutputWorkspace='diamond')