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
            ax = aLine[0]
            ay = aLine[1]
            az = aLine[2]
            aList = [ax,ay,az]
    return (V,aList)
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
"""
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

"""