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

def getEnergy(file): # this function returns E(sigma->0), not TOTEN
    """ parses an OUTCAR file and pulls out the energy of the system """
    f = open(file,'r')
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' in nextLine:
            f.readline()    # line of dashes
            f.readline()    # TOTEN line
            f.readline()    # blank line
            energyLine = f.readline().split()
            energy = float(energyLine[6])
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
            aLine = f.readline().split()
            ax = float(aLine[0])
            ay = float(aLine[1])
            az = float(aLine[2])
            aList = [ax,ay,az]
    return (V,aList)
    
def parseResults(directory):
    """ finds the energy of each simulation in a given results directory """
    dirList = []
    EList = []
    VList = []
    aLists = []
    for dir in os.listdir(directory):
        for file in os.listdir(directory + dir):        
                if fnmatch.fnmatch(file,'OUTCAR'):
                    EList += [getEnergy(directory+dir+'/'+file)]
                    sizeData = getCellSize(directory+dir+'/'+file)
                    VList += [sizeData[0]]
                    aLists += [sizeData[1]]
                    dirList += [dir]     
    return [dirList,EList,VList,aLists] 
            
#===========================================================================
# MAIN PROGRAM
#===========================================================================
direct = raw_input('HOME or WORK: ') # where all the results folders are located
if 'h' in direct or 'H' in direct:
    RESULTS = HOME
else:
    RESULTS = WORK
jName = raw_input('Job name: ')
data = parseResults(RESULTS+jName+'_results/')
dirList = data[0]
EList = data[1]
VList = data[2]
aLists = data[3]
print '\nDirectory names:\n'+str(dirList)
print '\nVolumes:\n'+str(VList)
print '\nLattice parameters:\n'+str(aLists)
perps = []
for a in aLists:
    perps += [a[1]]
print '\nPerpendicular parameters:\n'+str(perps)
print '\nEnergies:\n'+str(EList)

