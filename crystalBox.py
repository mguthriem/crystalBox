# some useful generation/packaging of crystallographic info using mantid 
import csv
import os
from mantid.simpleapi import *
# import matplotlib.pyplot as plt
import numpy as np
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter

class Box:

    '''class to hold list of peaks and their properties'''

    def __init__(self,cif):
 
        self.defaultCifFolder = '/SNS/SNAP/shared/cifLibrary'
        self.cifSpec = cif
        self.nickName = ''

        #locate and validate cif path using cifSpec.
        self.validCif = self.findCif()
        if not self.validCif:
            print('ERROR: validCif failed')

        self.dMin = 0.5
        self.dMax = 100

        #use mantid to derive crystallographic properties of peaks from cif
        if self.validCif:
            self.loadCif()


    def findCif(self):


        dirname = os.path.dirname(self.cifSpec)

        if self.cifSpec[-3:].lower() != 'cif':
            # look for a standard material: one that's listed by nickname in file cifIndex.csv 
            # stored in defaultCifFolder
            nickName=[]
            cifFilename=[]
            with open(f"{self.defaultCifFolder}/cifIndex.csv", mode = 'r') as file:
                csvFile = csv.reader(file)
                for line in csvFile:
                    nickName.append(line[0].lower()) #enforce lower case
                    cifFilename.append(line[1])
            try:
                foundIndex = nickName.index(self.cifSpec.lower())
                self.cifFilePath=f"{self.defaultCifFolder}/{cifFilename[foundIndex]}"
                self.nickName=self.cifSpec.lower()
                
            except:
                print(f"ERROR: cif specification {self.cifSpec} diamond failed")
                return False
            
            if os.path.isfile(self.cifFilePath):
                return True
            else:
                print(f"ERROR: tried and failed to open: {self.cifFilePath}")
                return False
            
        if (self.cifSpec[-3:].lower() == 'cif') and (dirname==''):

            #first see if file is in current working directory
            cwd =  os.path.dirname(__file__)
            cwdPath = f"{cwd}/{self.cifSpec}"
            if os.path.isfile(cwdPath):
                self.cifFilePath=cwdPath
                return True
            #next check defaultCifFolder
            
            libPath = f"{self.defaultCifFolder}/{self.cifSpec}"
            if os.path.isfile(libPath):
                self.cifFilePath = libPath
                return True
        
            print("ERROR: couldn\'t read cif file after checking these locations")
            print(cwdPath)
            print(libPath)
            return False

        if (self.cifSpec[-3:].lower() == 'cif') and (dirname != ''):

            if os.path.isfile(self.cifSpec):
                self.cifFilePath = self.cifSpec
                return True
            else:
                print("ERROR couldn\'t read this file: {self.cifSpec}")        
            
    def loadCif(self):
        CreateSampleWorkspace(OutputWorkspace='tmp')
        LoadCIF(Workspace='tmp',InputFile=self.cifFilePath)
        ws = mtd['tmp']
        crystal = ws.sample().getCrystalStructure()
        unitCell = crystal.getUnitCell()

        # define useful unit cell parameters
        self.a = unitCell.a()
        self.b = unitCell.b()
        self.c = unitCell.c()
        self.alpha = unitCell.alpha() 
        self.beta = unitCell.beta() 
        self.gamma = unitCell.gamma()
        self.volume = unitCell.volume()

        #keep a copy of lattice parameters in case they get scaled
        self.a_orig = unitCell.a()
        self.b_orig = unitCell.b()
        self.c_orig = unitCell.c()
        # some symmetry parameters

        self.HMSymbol = crystal.getSpaceGroup().getHMSymbol()
        
        #Generate reflections
        generator = ReflectionGenerator(crystal)
        # Create list of unique reflections between 0.7 and 3.0 Angstrom
        hkls = generator.getUniqueHKLsUsingFilter(self.dMin, self.dMax, ReflectionConditionFilter.StructureFactor)
        # Calculate d and F^2
        dValues = generator.getDValues(hkls)
        fSquared = generator.getFsSquared(hkls)
        self.pointGroup = crystal.getSpaceGroup().getPointGroup()

        # Make list of tuples and sort by d-values, descending, include point group for multiplicity.
        reflections = sorted([(hkl, d, fsq, len(self.pointGroup.getEquivalents(hkl))) for hkl, d, fsq in zip(hkls, dValues, fSquared)],
                                    key=lambda x: x[1] - x[0][0]*1e-6, reverse=True)

        # create individual lists of reflection properties with shared order and useful names

        self.nRef = len(reflections)
        self.hkl = []
        self.dspacing = []
        self.fsq = []
        self.mult = []
        self.estInt = []
        for i in range(self.nRef):
            self.hkl.append(reflections[i][0])
            self.dspacing.append(reflections[i][1])
            self.fsq.append(reflections[i][2])
            self.mult.append(reflections[i][3])
            Amp = reflections[i][2]*reflections[i][3]*reflections[i][1]**4 #Fsq times multiplicity * d**4
            self.estInt.append(Amp)

        DeleteWorkspace(Workspace='tmp')

    def summary(self):
        print(f"\nCIF file: {self.cifFilePath}")
        print(f"phase nickname: {self.nickName}")
        print(f"Space Group: {self.HMSymbol}")
        print(f"a: {self.a:.4f} Ang, b: {self.b:.4f} Ang, c: {self.c:.4f} Ang")
        print(f"alp: {self.alpha:.1f} deg, beta: {self.beta:.1f} deg, gam: {self.gamma:.1f} deg")
        print(f"{self.nRef} reflections calculated")
        print(f"First 10 reflections:")
        for ref in range(self.nRef):
            print(f"{self.hkl[ref]} {self.dspacing[ref]:4f} {self.mult[ref]} {self.estInt[ref]:.4f}")

    def tickWS(self,yVal):
        
        self.tickWSName = f"ticks: {self.nickName}"

        dataXArray = np.array(self.dspacing)
        dataYArray = np.ones_like(dataXArray)*yVal
        CreateWorkspace(OutputWorkspace=self.tickWSName,
                        DataX = dataXArray,
                        DataY = dataYArray,
                        UnitX = 'd-Spacing')
        
    def scaleLattice(self,scale):

        self.a = self.a_orig*scale
        self.b = self.b_orig*scale
        self.c = self.c_orig*scale

        return
    
    def resetLattice(self):

        self.a = self.a_orig
        self.b = self.b_orig
        self.c = self.c_orig

        return