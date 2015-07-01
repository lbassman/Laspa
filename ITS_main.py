#!/usr/bin/env python

# Jonas Kaufman jlkaufman@hmc.edu
# April 17, 2015
# Main script to run all calculations necessary to for finding
# the ideal tensile strength of a cubic crystal structure - on Stampede
# Building on work done by Josh Sanz

### TODO: clean up variable names


"""
First do a full relaxation, then move the OUTCAR and CONTCAR files to the
working directory as OUTCAR.E0 and CONTCAR.E0.

Add the Cell.py, POTCAR, KPOINTS, and INCAR.isif2 and INCAR.static files to the
working directory. (Also add the other python scripts, make them executable)

Make sure this script is executable.

Run this script to generate the necessary files and submit one strain step
relaxation job to the queue on Stampede.
"""

import subprocess as sp
from Cell import *
import os
import sys
import numpy as np

# number of cores for each calculation
NCORES = 16 
runLength = 100 ## to be set by user later on
# home and work directories
HOME = '/home1/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'
# Stampede allocation number
ALLOCATION = 'TG-DMR140093'
# email address for Slurm notifications
EMAIL = 'jlkaufman@hmc.edu'

def strainCell(cell, Epsilon):
    """
    applies a strain matrix Epsilon = [e1,e2,e3,e4,e5,e6] to the simulation
    cell
    """
    e1,e2,e3,e4,e5,e6 = Epsilon
    lvecs = cell.latticeVectors
    # modify lattice vectors
    for i in range(3):
        x,y,z = lvecs[i]
        lvecs[i] = [(1+e1)*x+(e6*y)/2+(e5*z)/2, \
                (e6*x)/2+(1+e2)*y+(e4*z)/2,(e5*x)/2+(e4*y)/2+(1+e3)*z]
    cell.setLatticeVectors(lvecs)
    # if in cartesian coordinates => modify atomic sites as well
    if cell.CorD[0] == 'C' or cell.CorD[0] == 'c':
        for i in range(cell.numberOfElements()):
            for j in range(cell.numberOfAtomsOfElement(i)):
                x,y,z = cell.sites[i][j].position
                cell.sites[i][j].move([(1+e1)*x+(e6*y)/2+(e5*z)/2,
                                       (e6*x)/2+(1+e2)*y+(e4*z)/2,
                                       (e5*x)/2+(e4*y)/2+(1+e3)*z])
    return cell

def genSubScript(jName,aList,runLength,NCORES):
    """
    create a submission script for Stampede's SLURM queueing system
    this version sends an email whenever one of the queued jobs starts
    """
    # uses integer division
    hrs = runLength/60
    mins = runLength%60
    string = ('#!/bin/bash\n' +
    '#SBATCH -J ' + jName +  '\n' +             # specify job name
    '#SBATCH -o ' + jName + '%j\n' +            # write output to this file
    '#SBATCH -n %d\n'%(NCORES*len(aList)) +     # request cores
    '#SBATCH -p normal\n' +                     # send to normal queue
    '#SBATCH -t %02d:%02d:00\n'%(hrs,mins) +      # set maximum wall (clock) time
    '#SBATCH --mail-user=' + EMAIL +'\n' +      # set email
    '#SBATCH --mail-type=all\n' +               # send all emails
    '#SBATCH -A ' + ALLOCATION + '\n' +         # specifies project
    'module load vasp\n')                       # load vasp module
    for i in range(len(aList)):
        # change to work directory, run vasp
        string += 'cd '+WORK+'%s_%.5f\n'%(jName,aList[i])
        string += "ibrun -o %d -n %d vasp_std > vasp_output.out &\n"%(NCORES*i,NCORES)
    string += 'wait\ncd '+HOME+'\nmkdir %s_results\n'%(jName)
    for i in range(len(aList)):
        # move directories to results directory
        string += 'cd '+WORK+'\nmv %s_%.5f '%(jName,aList[i])
        string += HOME+'%s_results/\n'%jName
    f = open(jName + '_submit','w')
    f.write(string)
    f.close()

def runITS(jName,aList,runLength,NCORES):
    """
    generate subdirectory, copies files, and submits job to
    relax a strained cell
    """
    # copy files to subdirectory, move subdirectory to WORK
    for a in aList:
        # create submission script
        genSubScript(jName,aList,runLength,NCORES)
        # apply perpendicular strains
        cell = Cell().loadFromPOSCAR('POSCAR.EN')
        strainCell(cell, [0,a,a,0,0,0])
        cell.sendToPOSCAR()
        # copy files to subdirectory, move subdirectory to WORK
        sp.call(['mkdir','%s_%.5f'%(jName,a)])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+[jName+'_submit']+\
                ['%s_%.5f'%(jName,a)])
        sp.call('cp -r %s_%.5f '%(jName,a)+WORK,shell=True)
        sp.call('chmod u+x %s_submit'%jName,shell=True)
    # run submission script
    sp.call(['sbatch','%s_submit'%jName])

#===========================================================================
# MAIN PROGRAM
#===========================================================================
NCORES = float(raw_input('Number of cores per simulation: '))
emax = float(raw_input('Maximum strain (%): '))/100.0     # e.g. 15%
inc = float(raw_input('Strain step size (%): '))/100.0    # e.g. 1%
NLAT = int(raw_input('Number of perpendicular lattice vector strains to run: '))
job = raw_input('Job name: ')
valid = 0
while not valid:
    try:
        runLength = int(raw_input("Maximum run time in minutes: "))
        valid = 1
    except:
        print "Run time must be a whole number"
eList = np.arange(inc,emax+inc,inc)
for eN in eList:
    aList = np.linspace(0.0, -eN, NLAT)
    # apply strain
    cell = Cell().loadFromPOSCAR('CONTCAR.E0') # starts from provided POSCAR
    strainCell(cell,[eN,0,0,0,0,0])
    cell.sendToPOSCAR('POSCAR.EN')    
    jName = job+'_%.3f'%eN
    # run VASP job
    print 'Running %.3f strain'%eN
    runITS(jName,aList,runLength,NCORES)