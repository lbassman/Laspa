#!/usr/bin/env python

# Jonas Kaufman jlkaufman@hmc.edu
# June 29, 2015
# Script to parse VASP output directories for energy values

import os
import fnmatch
import numpy as np

# home and work directories
HOME = '/home1/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'

def getEnergy(file): # change to parse OSZICAR?
    """ parses an OUTCAR file and pulls out the free energy of the system """
    f = open(file,'r')
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' in nextLine:
            f.readline()    # line of dashes
            energyLine = f.readline().split()
            energy = float(energyLine[4])
    return energy
    
def getCellSize(file):
    """ returns the volume of a cell parsed from OUTCAR """
    f = open(file,'r')
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'VOLUME and BASIS-vectors are now :' in nextLine:
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            volumeLine = f.readline().split()
            V = float(volumeLine[4])
            for i in range(6):
                f.readline()    # text
            a = float(f.readline().split()[1]) # perpendicular lattice vectors
            # the lattice vector length in loading direction is [0] ## change this?
            break
    return (V,a)
    
def parseResults(directory):
    """ finds the energy of each simulation in a given results directory """
    dirList = []
    EList = []
    VList = []
    for dir in os.listdir(directory): # 0.XXXXX_main or 0.XXXXX_c*_0.XXXXX
        for file in os.listdir(directory + dir):        
                if fnmatch.fnmatch(file,'OUTCAR'):
                    EList += [getEnergy(directory+dir+'/'+file)]
                    sizeData = getCellSize(directory+dir+'/'+file)
                    VList += [sizeData[0]]
                    # aList += [sizeData[1]]
                    dirList += [dir]     
    return [dirList,EList,VList] 
            
#===========================================================================
# MAIN PROGRAM
#===========================================================================
RESULTS = HOME # where all the results folders are located
cname = raw_input('Job name: ')
data = parseResults(RESULTS+cname+'_results/')
dirList = data[0]
EList = data[1]
VList = data[2]
print '\nDirectory names:\n'+str(dirList)
print '\nVolumes:\n'+str(VList)
print '\nEnergies:\n'+str(EList)
