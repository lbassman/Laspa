#!/usr/bin/env python

# Jonas Kaufman jlkaufman@hmc.edu
# June 10, 2015
# Helper script to parse VASP output files for data necessary for
# finding the ideal tensile strength of a cubic crystal structure - on Stampede
# Building on work done by Josh Sanz

import os
import fnmatch
import numpy as np

# home and work directories
HOME = '/home1/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'

def readState():
    # open state file and get information
    #decimals, not percents
    state = open('STATE','r+')
    line1 = state.readline()
    if line1 == '':
        eN = 0.0
        inc = INC
        emax = EMAX
    else:
        eN = float(line1.split()[1])
        inc = float(state.readline().split()[1])
        emax = float(state.readline().split()[1])
    state.close()
    return [eN, inc, emax]

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
    
def getInitialConditions():
    """ collects the equilibrium energy, volume, and lattice constant from
        the OUTCAR.E0 file """
    E0 = getEnergy('OUTCAR.E0')
    V0,a0 = getCellSize('OUTCAR.E0')
    return (E0,V0,a0)

def parseResults(directory):
    """ finds the energy of each simulation in a given results directory """
    EList = []
    VList = []
    aList = []
    tList = []
    for dir in os.listdir(directory): # 0.XXXXX_main or 0.XXXXX_c*_0.XXXXX
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
# read in state file to get strain steps
state = readState()
inc = state[1]
emax = state[2]
eList = []
eN = 0.0
while abs(eN) < (abs(emax)+inc/2.0): # fix this?
    eList += [eN]
    eN += inc
print eList

RESULTS = HOME # where all the results folders are located
for eN in eList:
    # get data from main strain steps
    print eN
    print 'main strain results:'
    if eN == 0.0:
        # get initial conditions of completely unstrained cell
        #print getInitialConditions()+[0.0]
        print 'NEED to add OUTCAR.E0'
    else:
        print parseResults(RESULTS+'%.5f_main_results/'%eN)
    # get data from additional strain runs
    print 'cprime results:'
    print parseResults(RESULTS+'%.5f_cprime_results/'%eN)
    print 'c44 results:'
    print parseResults(RESULTS+'%.5f_c44_results/'%eN)
    print 'c66 results:'
    print parseResults(RESULTS+'%.5f_c66_results/'%eN)

#strains += [(a-a0)/a0] ?