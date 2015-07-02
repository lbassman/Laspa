#!/usr/bin/env python

# Jonas Kaufman jlkaufman@hmc.edu
# June 29, 2015
# Script to parse VASP output directories
#### changing this to be more of a "dashboard" for analyzing VASP runs

# Run 'module load python' before using on Stampede
# Make script executable 'chmod +x _____.py'

import os
import fnmatch
import numpy as np
import matplotlib as plt

# home and work directories
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'

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

def getEnergies(file): # this function returns E(sigma->0), not TOTEN
    """ parses an OUTCAR file and pulls out the energy of the system
    after each ionic step """
    energies = []
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
            energies += [float(energyLine[6])]
    return energies
    
def getCellSize(file):
    """ returns the volume and lattice parameters """
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

def getSizes(file):
    """ returns the volume and lattice parameters
    after each ionic step
    """
    f = open(file,'r')
    volumes = []
    lats = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'VOLUME and BASIS-vectors are now :' in nextLine:
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            volumeLine = f.readline().split()
            volumes += [float(volumeLine[4])]
            for i in range(6):
                f.readline()    # text
            aLine = f.readline().split()
            ax = float(aLine[0])
            ay = float(aLine[1])
            az = float(aLine[2])
            lats += [ax,ay,az]
    return (volumes, lats)    
    
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
    return (dirList,EList,VList,aLists) 
            
#===========================================================================
# MAIN PROGRAM
#===========================================================================
"""
RESULTS = HOME # where all the results folders are located # change to be input
jName = raw_input('Job name: ')
data = parseResults(RESULTS+jName+'_results/')
dirList = data[0]
EList = data[1]
VList = data[2]
aLists = data[3]
print '\nDirectory names:\n'+str(dirList)
print '\nVolumes:\n'+str(VList)
print '\nLattice parameters:\n'+str(aLists)
print '\nEnergies:\n'+str(EList)
"""
plt.plot([1,2,3,4])
plt.ylabel('some numbers')
plt.show()
