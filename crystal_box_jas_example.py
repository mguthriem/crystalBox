# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import importlib

from mantid.geometry import PointGroupFactory
from mantid.geometry import SpaceGroupFactory

sys.path.append("/SNS/users/xkq/Documents/SNAP_SXtal_software_dev/crystalBox/")
import crystalBox as crys
importlib.reload(crys)

#1. Create a "crystalBox" using a full path to a cif file
silicon = crys.Box('/SNS/users/xkq/Documents/SNAP_SXtal_software_dev/cifs/Si_2001_Toebbens.cif')

#2. Display summary of crystal information using

print("A summary of the silicon crystal info:")
print(silicon.summary())

#3. See if the function I made works
print ("This is a test of if the function I made works:")
print(silicon.jasmineFunction())

print("This is a test of my 'getEquivalents(self,hkl)' function:")
print(silicon.getEquivalents([2,2,0]))

print("This is a test of new method 'makeCrystal':")
crystal = silicon.makeCrystal()
print(crystal)

print("I think I need to run malcolm's 'proccessCrystal' method first:")
print(silicon.processCrystal(crystal))

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