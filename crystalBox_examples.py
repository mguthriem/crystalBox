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

nickel.reset()

print(f"reset nickel lattice param: {nickel.a:.4f} Ang, vol = {nickel.volume:.4f} Ang^3")
print(nickel.dSpacing)

#6. You can create a workspace with ticks using

print(f"\nNow creating workspace \'ticks: diamond\' with diamond peak positions")
diamond.tickWS(-0.01) #the value in brackets is the height of the tick marks


#you can then drag diffraction data onto that plot

LoadNexus(Filename='/SNS/SNAP/shared/Malcolm/code/crystalBox/data/diamond.nxs',
    OutputWorkspace='diamond')
    
#and then plot is using

diamond.plot("diamond")
    
#Jasmine's examples:

#2. Display summary of crystal information using

silicon = crys.Box('silicon')

print("A summary of the silicon crystal info:")
print(silicon.summary())

print("This is a test of my 'getEquivalents(self,hkl)' function:")
equivs = silicon.getEquivalents([2,2,0])
print(equivs)


print("This is a test of the function 'cartesianHKL':")
print(silicon.cartesianHKL([4.64539, -1.34246, 3.24761]))

print("This is a test of the function 'cartesianHKL':")
print(silicon.cartesianHKL([4.66322, -3.16803, -1.46518]))

print("This is a test of the function 'getAngle':")
angle_calc = silicon.getAngle([2, 5, 1], [-7, 1, -4])
print(angle_calc)

print("Now retreiving angle between corresponding q-coordinates:")
angle_obs = silicon.getAngle([4.64539, -1.34246, 3.24761], [-6.83164, 1.4796, 6.20429])
print(angle_obs)

print("This is a test of the function 'getAngle' where hkls are not rounded:")
angle_calc = silicon.getAngle([1.96187, 4.98722, 0.91789], [-6.95302, 0.958652, -3.97906])
print(angle_calc)

print(silicon.defineSpaceGroup())