# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import importlib

sys.path.append("/SNS/SNAP/shared/Malcolm/code/crystalBox")
import crystalBox as crys
importlib.reload(crys)

#1. Create a "crystalBox":
crystal = crys.Box('kdp')

print("\nTEST 1: create original crystal Box and display \n")
crystal.summary()

#2. Scale Lattice and check this works
print("\nTEST 2: Scaling lattice crystal \n")
crystal.scaleLattice(0.99)
#check this worked
if crystal.isModified:
    print("crystal structure has been modified")
    crystal.summary()
else:
    print("Crystal is unmmodified")
    crystal.summary

#3. Resetting crystal and check this works
print("\nTEST 3:resetting crystal \n")
crystal.reset()
#check again
if crystal.isModified:
    print("crystal structure has been modified")
    crystal.summary()
else:
    print("Crystal is unmmodified")
    crystal.summary

#4. Modify Uiso and check this works

print("\nTEST 4: convert original uiso to \"standardised value\" \n")
crystal.applyStandardUiso()
#check again
if crystal.isModified:
    print("crystal structure has been modified")
    crystal.summary()
else:
    print("Crystal is unmmodified")
    crystal.summary

#3. Again Reset crystal and check this works
print("\nTEST 5: resetting crystal \n")
crystal.reset()
#check again
if crystal.isModified:
    print("crystal structure has been modified")
    crystal.summary()
else:
    print("Crystal is unmmodified")
    crystal.summary