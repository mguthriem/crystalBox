import csv
import os
from mantid.simpleapi import *
import numpy as np
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter
import matplotlib.pyplot as plt
from mantid.plots.utility import MantidAxType
from mantid.api import AnalysisDataService as ADS
from mantid.kernel import Atom

from mantid.geometry import PointGroupFactory
from mantid.geometry import SpaceGroupFactory
#this is a test by Jasmine K. Hinton# 

defaultCifFolder = '/SNS/SNAP/shared/cifLibrary'

class parseMantidCrystal():

    '''an interface to mantid's CrystalStructure objects
    
    The Mantid CrystalStructure class 
    https://github.com/mantidproject/mantid/blob/main/Framework/Geometry/inc/MantidGeometry/Crystal/CrystalStructure.h
    
    Are instantiated using three key piece of crystallographic information:

    1. The lattice paramaters, specified as a string containing 3 or 6 floating point numbers (if 3 are used, the final 3 are taken to be 90.0 90.0 90.0)
    2. The space group, also specified as a string
    3. The atomic basis, specified as a single string where each atom has space separated
        Type x y z occupancy uiso, and each atom is separated by a semi colon.

    Meanwhile these same quantities are output by dedicated methods in the following way:

    1. getUnitCell() return a class that contains the individual lattice parameters as individual
        individual attributes.
    2. getSpaceGroup().getHMSymbol() is a string
    3. getScatterers() returns a list of string, one string per atom
    
    The goal of localCrystalStructure is to provide a more convenient holder for the crystallographic 
    properties that allows them to be easily modified.

    It will be instantiated by a mantid CrystalStructure object and has the following attributes:

    1. unitCellList = a list of floats corresponding to a,b,c,alpha,beta,gamma
    2. spaceGroupString = is identical to getSpaceGroup().getHMSymbol()
    3. cellContentsList = a list of lists, with each atom being represented by a list of its properties

    It will also have methods to convert its own attributes into a mantid CrystalStructure.

    '''

    def __init__(self,mantidCrystalStructure):

        #useful attributes stored as lists...

        #unitCellList
        unitCellObject = mantidCrystalStructure.getUnitCell()
        self.unitCellList = [unitCellObject.a(),
                   unitCellObject.b(),
                   unitCellObject.c(),
                   unitCellObject.alpha(),
                   unitCellObject.beta(),
                   unitCellObject.gamma()]
        
        #cellContentsList
        cellContentsString = mantidCrystalStructure.getScatterers()
        cellContentsList = [atm.split(' ') for atm in list(cellContentsString)]
        self.cellContentsList = [[val if val.isalpha() else float(val) \
                       for val in scatterer] for scatterer in cellContentsList]

        #spaceGroupString
        self.spaceGroupString = mantidCrystalStructure.getSpaceGroup().getHMSymbol()
        
    def makeMantidCrystal(self):

        #returns a mantid CrystalStructure 
        latticeParamsString = ''
        for param in self.unitCellList:
            latticeParamsString += f"{param} "

        cellContentsString = "; ".join([f'{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]} {atom[5]}' for atom in self.cellContentsList])
        
        mantidCrystalStructure = CrystalStructure(latticeParamsString,
                                         self.spaceGroupString,
                                         cellContentsString)

        return mantidCrystalStructure

class Box():

    '''class to hold list of peaks and their properties'''

    def __init__(self,cif):
 
        self.defaultCifFolder = defaultCifFolder
        self.cifSpec = cif
        self.nickName = ''

        #locate and validate cif path using cifSpec.
        self.validCif = self.findCif()
        if not self.validCif:
            print('ERROR: validCif failed')

        self.dMin = 0.5
        self.dMax = 100
 
        self.modifiedLattice = False # if false allows scaling of lattice
        if self.validCif:
            self.loadCif()           

        self.tickWSExists = False

        #plotting attributes
        self.markerHeight = 10
        self.markeredgecolor='red'


    def findCif(self):

        dirname = os.path.dirname(self.cifSpec)

        if self.cifSpec[-3:].lower() != 'cif':
            # look for a standard material: one that's listed by nickname in file cifIndex.csv 
            # stored in defaultCifFolder
            nickName=[]
            cifFilename=[]
            with open(f"{self.defaultCifFolder}/nickNames.csv", mode = 'r') as file:
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


        #Load cif into mantid workspace and get useful mantid object from this.
        CreateSampleWorkspace(OutputWorkspace='tmp')
        LoadCIF(Workspace='tmp',InputFile=self.cifFilePath)
        ws = mtd['tmp']
        
        #at point of loading make a copy of input crystal
        inputCrystalStructure =ws.sample().getCrystalStructure()
        inputCrystalList = parseMantidCrystal(inputCrystalStructure)

        #make a copy of input crystal as a mantid Crystal Structure
        self._originalCrystal = inputCrystalList.makeMantidCrystal()

        self.processCrystal(inputCrystalStructure)

        DeleteWorkspace(Workspace='tmp')

    def sameCrystal(self,crystalStructure1,crystalStructure2):

        #compares two mantid CrystalStructure objects and returns True if they have identical unitCells, spaceGroups and atomLists

        #first convert mantid CrystalStructure to crystalList
        crystal1 = parseMantidCrystal(crystalStructure1)
        crystal2 = parseMantidCrystal(crystalStructure2)

        #this was far harder than I though. Here I'm assuming the atoms never change their @order which seems not guaranteed. TODO: worry about this.
        scattererCondition = []
        for atomIndex,atom in enumerate(crystal1.cellContentsList):
            atomSet1 = set(atom)
            atomSet2 = set(crystal2.cellContentsList[atomIndex])
            scattererCondition.append(atomSet1 == atomSet2)
        scattererCondition = all(scattererCondition)

        spaceGroupCondition = crystal1.spaceGroupString == crystal2.spaceGroupString

        latticeSet1 = set(crystal1.unitCellList)
        latticeSet2 = set(crystal2.unitCellList)
        latticeCondition = latticeSet1 == latticeSet2

        if scattererCondition and spaceGroupCondition and latticeCondition:
            return True
        else:
            return False
    
    def processCrystal(self,crystalStructure):
        
        #Accepts a mantid Crystal Structure object and extracts useful attributes from it

        # first check if crystal has been modified from original as read from cif
        self.isModified = not self.sameCrystal(crystalStructure,
                                                      self._originalCrystal)

        #make lists of useful crystal parameters by parsing crystalStructure
        self.crystalList = parseMantidCrystal(crystalStructure)
        self.unitCellList = self.crystalList.unitCellList
        self.cellContentsList =  self.crystalList.cellContentsList
        self.nAtoms = len(self.cellContentsList)

        # some symmetry parameters
        self.HMSymbol = crystalStructure.getSpaceGroup().getHMSymbol()
        self.pointGroup = crystalStructure.getSpaceGroup().getPointGroup()
        self.crystalSystem = str(self.pointGroup.getCrystalSystem())
        
        #define useful crystallographic attributes
        unitCell = crystalStructure.getUnitCell()  
        self.a = unitCell.a()
        self.b = unitCell.b()
        self.c = unitCell.c()
        self.alpha = unitCell.alpha() 
        self.beta = unitCell.beta() 
        self.gamma = unitCell.gamma()
        self.volume = unitCell.volume()
        self.astar = unitCell.astar()
        self.bstar = unitCell.bstar()
        self.cstar = unitCell.cstar()
        self.alphastar = unitCell.alphastar()
        self.betastar = unitCell.betastar()
        self.gammastar = unitCell.gammastar()

        #JKH added to make use in cartesian hkl method
        self.BMatrix = unitCell.getB() 

        #Generate reflections
        generator = ReflectionGenerator(crystalStructure)
        
        # Create list of unique reflections between 0.7 and 3.0 Angstrom
        hkls = generator.getUniqueHKLsUsingFilter(self.dMin, self.dMax, ReflectionConditionFilter.StructureFactor)
        # Calculate d and F^2
        dValues = generator.getDValues(hkls)
        fSquared = generator.getFsSquared(hkls)
        
        # Make list of tuples and sort by d-values, descending, include point group for multiplicity.
        reflections = sorted([(hkl, d, fsq, len(self.pointGroup.getEquivalents(hkl))) for hkl, d, fsq in zip(hkls, dValues, fSquared)],
                                    key=lambda x: x[1] - x[0][0]*1e-6, reverse=True)

        
        # create individual lists of reflection properties with shared order and useful names
        self.nRef = len(reflections)
        self.hkl = []
        self.dSpacing = []
        self.fSq = []
        self.mult = []
        self.estInt = []
        self.totalIntensityInRange = 0
        for i in range(self.nRef):
            self.hkl.append(reflections[i][0])
            self.dSpacing.append(reflections[i][1])
            self.fSq.append(reflections[i][2])
            self.mult.append(reflections[i][3])
            Amp = reflections[i][2]*reflections[i][3]*reflections[i][1]**4 #Fsq times multiplicity * d**4
            self.estInt.append(Amp)
            self.totalIntensityInRange += Amp/self.volume

        return
        
    def summary(self):
        print("Crystal summary:")
        print(f"\nCIF file: {self.cifFilePath}")
        print(f"phase nickname: {self.nickName}")
        print(f"Space Group: {self.HMSymbol}")
        print(f"a: {self.a:.4f} Ang, b: {self.b:.4f} Ang, c: {self.c:.4f} Ang")
        print(f"alp: {self.alpha:.1f} deg, beta: {self.beta:.1f} deg, gam: {self.gamma:.1f} deg")
        print(f"{self.nRef} reflections calculated")
        print(f"Atoms:")
        for atom in self.cellContentsList:
            print(atom)
        print("Crystal has been modified: ",self.isModified)
        print(f"First 3 reflections:")
        for ref in range(3):
            print(f"{self.hkl[ref]} {self.dSpacing[ref]:4f} {self.mult[ref]} {self.fSq[ref]:.4f}")

    def tickWS(self,yVal):

        self.tickWSExists = True
        self.tickWSName = f"ticks: {self.nickName}"
        self.tickWSyVal = yVal

        dataXArray = np.array(self.dSpacing)
        dataYArray = np.ones_like(dataXArray)*self.tickWSyVal
        CreateWorkspace(OutputWorkspace=self.tickWSName,
                        DataX = dataXArray,
                        DataY = dataYArray,
                        UnitX = 'd-Spacing')
        
    def scaleLattice(self,scaleFactor):

        print("updating lattice parameters")

        #scaleLattice: crude multiplication of original input a,b,c
        self.latticeScaleFactor = scaleFactor

        originalCrystalList = parseMantidCrystal(self._originalCrystal)
        originalUnitCellList = originalCrystalList.unitCellList

        #modify current crystalList with updated lattice parameters
        self.crystalList.unitCellList = [scaleFactor*val for val in originalUnitCellList]

        #build Crystal Structure from this
        crystalMod=self.crystalList.makeMantidCrystal()
        #and process this to update current values
        self.processCrystal(crystalMod)

        #if workspace exists, need to update it with scaled d-spacings
        if self.tickWSExists:
            self.tickWS(self.tickWSyVal)

        return
    
    def calcFOM(self):

        #FOM is a generalised figure of merit assessing how crystallographically #challenging a material is

        self.dMin = 0.5
        self.dMax = 5.0
        self.applyStandardUiso()
        FOM = self.totalIntensityInRange/self.nRef
        self.FOM1 = FOM/548.881 #normalised to value for diamond

        # print(f"FOM is: {self.FOM1}")

        #a second FOM assessess total absorption and incoherent cross sections.
        AIFOM = self.calcIncAbsFactor()
        self.FOM2 = AIFOM/0.0007931 #normalised to value for diamond
        # print(f"Abs Inc FOM is: {self.FOM2}")
        return self.FOM1,self.FOM2
        

    def applyStandardUiso(self):

        # print("updating Uiso values...")


        for atom in self.crystalList.cellContentsList:
            Mass = Atom(atom[0]).mass
            #factor of 30 empirically chosen to make the F2 numbers less vanishingly small
            atom[5] = 1/(30*np.sqrt(Mass))

        crystalMod = self.crystalList.makeMantidCrystal()
        self.processCrystal(crystalMod)
        
        return
    
    def calcIncAbsFactor(self):
        #Calculates a numerical factor to assist in determining how 'well' as sample scatters
        #A large returned number from this calculator will scatter less well than a lower number structure.
        contents=self.cellContentsList
        sg = SpaceGroupFactory.createSpaceGroup(self.HMSymbol)
        factor = 0
        for i in range(0,len(contents)):
            xCoord = contents[i][1]
            yCoord = contents[i][2]
            zCoord = contents[i][3]
            siteOcc = contents[i][4]
            position = [xCoord, yCoord, zCoord]
            siteMult = len(sg.getEquivalentPositions(position))
            atom = Atom(contents[i][0])
            incXS = atom.neutron()['inc_scatt_xs']
            absXS = atom.neutron()['abs_xs']
            #The factor is the sum of the site multiplicity multiplied by the occupancy and (inc+abs) factor for every atom in the cell.
            factor += siteMult * siteOcc * (incXS + absXS)
        #The factor is scaled by the cell volume. A larger cell with the same number of atoms has a lower density of scatterers (absorption or inc).
        factor = factor / self.volume
        return factor   
        
    def reset(self):

        #returns entire crystal to original as read from cif
        self.processCrystal(self._originalCrystal)
        return

    
    def plot(self,workspaceToPlot):
                        
        """workspaceToPlot is the an existing workspace that will be plotted along wit the ticks 
        """

        ticks = ADS.retrieve(self.tickWSName)
        if workspaceToPlot != '':
            dataToPlot = ADS.retrieve(workspaceToPlot)
            nHst = mtd[workspaceToPlot].getNumberHistograms()
            if nHst >= 7:
                nHst=6
                print("WARNING: can only plot first 6 spectra!")


        fig, axes = plt.subplots(edgecolor='#ffffff', num='ticks plot', subplot_kw={'projection': 'mantid'})
        lineColours = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']

        axes.plot(ticks, color='#1f77b4', 
                  label=f"ticks: {self.nickName}", 
                  linestyle='None', 
                  marker='|', 
                  markersize=self.markerHeight,
                  markeredgecolor=self.markeredgecolor, 
                  wkspIndex=0)

        if workspaceToPlot != '':
            for spectrum in range(nHst):
                axes.plot(dataToPlot,wkspIndex=spectrum,color=lineColours[spectrum],label=f"spec{spectrum}")

        axes.tick_params(axis='x', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
        axes.tick_params(axis='y', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
        axes.set_title(f"ticks: {self.nickName}")
        axes.set_xlabel('d-Spacing ($d-Spacing$)')
        axes.set_ylabel('($d-Spacing$)$^{-1}$')
        # axes.set_xlim([0.4448, 2.0285])
        # axes.set_ylim([0.14175, 0.15825])
        legend = axes.legend(fontsize=8.0).set_draggable(True).legend

        plt.show()

    def dLimits(self,dMin,dMax):
        self.dMin=dMin
        self.dMax=dMax
        self.loadCif()

    
    #########################################################################
    #                           getEquivalents                              #
    #   getEquivalents is intended to take a reflection, pull the space     #
    #   group from the cif, create a point group, and then provide all the  #
    #   equivalent reflections of that given hkl.                           #
    #########################################################################

    def getEquivalents(self, hkl):

        equivalents = self.pointGroup.getEquivalents(hkl)

        return equivalents      
       
    #########################################################################
    #                           cartesianHKL                                #
    #   the purpose of cartesianHKL is to take a vector, in this case,      #
    #   likely a set of 3 q-coordinates, and if they're not cartesian, then #
    #   multiply by the b-matrix (obtained in process crystal method) and   #
    #   the resulting vector should then be cartesian. It now actually      #
    #   verifies the crystal system has a cartesian basis and proceeds with #
    #   the calculation if so. Otherwise it says you you can skip this step #                                                                   #
    #########################################################################
    
    def cartesianHKL(self,h):

        crystalSystem = str(self.pointGroup.getCrystalSystem())

        if crystalSystem == "Cubic" or crystalSystem == "Orthorhombic" or crystalSystem == "Tetragonal":
            print("Crystal system is ", crystalSystem)
            print("Crystal system has cartesian basis and does not need to be converted via cartesianHKL.",
                  "You may proceed to calculating angle if you wish.")
        else:
            print("Crystal system is ", crystalSystem)
            print("pg.getCrystalSystem() outputs a variable of type:", type(crystalSystem))
            print("Crystal system does not have cartesian basis and will now be converted via caresianHKL")
            
            v = np.matmul(self.BMatrix,h)
            return v
    
    #########################################################################
    #                           getAngle                                    #
    #   getAngle is intended to take two vectors and calculate the angle    #
    #   between them in degrees. A boolean variable called degrees is       #
    #   assigned to help the user know that the result will be in degrees   #
    #   and not radians. It is assumed that the vectors are on an           #
    #   an orthonormal cartesian basis. If they are not, it is recommended  #
    #   to use the method cartesian hkl first.                              #
    #########################################################################

    def getAngle(self,vector1,vector2, degrees=True):
        
        print("Ensure your vectors are in a cartesian basis set before using this method,",
              "otherwise your result will not be meaningful.",
              "Use method, 'cartesianHKL' to verify this.")
        
        # ensure vectors use cartesian metric

        if self.crystalSystem == "Cubic" or self.crystalSystem == "Orthorhombic" or self.crystalSystem == "Tetragonal":
            self.isCartesian = True
        else:
            self.isCartesianm = False
            vector1 = np.matmul(self.BMatrix,vector1)
            vector2 = np.matmul(self.BMatrix,vector2)
            print("NOTICE: vectors converted to cartesian")
        
        angleBetweenVectors = np.arccos(np.dot(vector1,vector2)/(np.linalg.norm(vector1)*
            np.linalg.norm(vector2)))
        
        if degrees:
            angleBetweenVectors = np.degrees(angleBetweenVectors)
        
        return angleBetweenVectors 
        
def showNicknames():
    print('available nicknames are:')
    with open(f"{defaultCifFolder}/nickNames.csv", mode = 'r') as file:
        csvFile = csv.reader(file)
        for line in csvFile:
            print(line[0].lower())
