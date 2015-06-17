# Jonas Kaufman jlkaufman@hmc.edu
# March 7, 2015
# Script to parse VASP output files for energy values
# for a set of lattice parameters for cubic structures
# Based on a script by Josh Sanz

import os
import fnmatch

# home and work directories
HOME = '/home1/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'

def getEnergy(file):
    """ parses an OUTCAR file and pulls out the free energy of the system """
    f = open(file,'r')
    while True:
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' in f.readline():
            f.readline()    # line of dashes
            energyLine = f.readline().split()
            energy = float(energyLine[4])
            break
    return energy
    
def getCellSize(file):
    """ returns the volume of a cell parsed from OUTCAR """
    f = open(file,'r')
    while True:
        if 'VOLUME and BASIS-vectors are now :' in f.readline():
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            volumeLine = f.readline().split()
            V = float(volumeLine[4])
            for i in range(6):
                f.readline()    # text
            a = float(f.readline().split()[1]) # perpendicular lattice vectors
            # the lattice vector length in loading direction is [0]
            break
    return (V,a)
    ## check that this is really reading the perpendicular size

def parseResults(directory):
    """ finds the energy of each simulation in a given results directory """
    EList = []
    VList = []
    aList = []
    tList = []
    for dir in os.listdir(directory): # 0.XXXXX_main or 0.XXXXX_c*_0.XXXXX for example
        for file in os.listdir(directory + dir):        
                if fnmatch.fnmatch(file,'OUTCAR'):
                    EList += [getEnergy(directory+dir+'/'+file)]
                    sizeData = getCellSize(directory+dir+'/'+file)
                    VList += [sizeData[0]]
                    aList += [sizeData[1]]
        dirInfo = dir.split('_')
        # eN = dirInfo[0] # main strain amount
        if len(dirInfo) == 2:
            tList += [0.0] # additional strain amount
        else:
            tList += [float(dirInfo[2])]          
    return [EList,VList,aList,tList] 


#===========================================================================
# MAIN PROGRAM
#===========================================================================
RESULTS = HOME # where all the results folders are located

print 'GS results:'
data = parseResults(RESULTS+'GS_results/')
EList = data[0]
VList = data[1]
print '\nVolumes:\n'+str(VList)
print '\nEnergies:\n'+str(EList)
print 

   