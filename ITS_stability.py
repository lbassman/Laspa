#!/usr/bin/env python

# Jonas Kaufman jlkaufman@hmc.edu
# June 10, 2015
# Helper script to apply additional strains to test for stability and
# run VASP calculations in order to find the ideal tensile strength
# of a cubic crystal structure - on Stampede
# Building on work done by Josh Sanz

"""
First do a full relaxation, then move the OUTCAR and CONTCAR files to the
working directory as OUTCAR.E0 and CONTCAR.E0.

Generate OUTCAR and CONTCAR files for each strain step using the main script
and make sure they are in the same directory as CONTCAR.E0.XXXXX (or .E-0.XXXXX)

Add the Cell.py, POTCAR, KPOINTS, and INCAR.isif2 and INCAR.static files to the
working directory.

Run this script to generate the necessary files submit stability checking jobs
to the queue on on Stampede
"""

import subprocess as sp
from Cell import *
import os

# number of cores for each calculation
NCORES = 16
runLength = 5 ## to be set by user later on

# home and work directories
HOME = '/home1/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'

# email address for Slurm notifications
EMAIL = 'jlkaufman@hmc.edu'

# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

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

def genSubScript(cname,tList,runLength,NCORES):
    """
    create a submission script for Stampede's SLURM queueing system
    this version sends an email whenever one of the queued jobs starts and ends
    """
    # uses integer division
    hrs = runLength/60
    mins = runLength%60
    string = ('#!/bin/bash\n' +
    '#SBATCH -J ' + cname +  '\n' +             # specify job name
    '#SBATCH -o ' + cname + '%j\n' +            # write output to this file
    '#SBATCH -n %d\n'%(NCORES*len(tList)) +     # request 64 cores
    '#SBATCH -p normal\n' +                     # send to normal queue
    '#SBATCH -t %.2d:%d:00\n'%(hrs,mins) +      # set maximum wall (clock) time
    '#SBATCH --mail-user=' + EMAIL +'\n' +      # set email
    '#SBATCH --mail-type=begin\n' +             # send email when job starts
    '#SBATCH --mail-type=end\n' +               # send email when job ends
    '#SBATCH -A ' + ALLOCATION + '\n' +         # specifies project
    'module load vasp\n')                       # load vasp module
    for i in range(len(tList)):
        # change to working directory, run vasp
        string += 'cd '+WORK+cname+'_%.5f\n'%tList[i]
        string += 'ibrun -o %d -n %d vasp_std > vasp_output.out &\n'%(NCORES*i,NCORES)
    # wait for all vasp runs to finish, move all files to home directory
    string += 'wait\ncd '+HOME+'\nmkdir %s_results\n'%(cname)
    for i in range(len(tList)):
        # move directories to results directory
        string += 'cd '+WORK+'\nmv '+cname+'_%.5f '%tList[i]+HOME+'%s_results/\n'%cname
    f = open(cname + '_submit','w')
    f.write(string)
    f.close()

def getCxx(cname,eN,tList,runType,runLength,NCORES):
    """
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK directory and runs submission scripts ???
    """
    if eN == 0.0:
        CONTCAR = 'CONTCAR.E0'
    else:
        CONTCAR = HOME+'%.5f_main_results/'%eN+'%.5f_main/CONTCAR'%eN
    for t in tList:
        # generate lattice w/ applied strain
        # start from equilibrium lattice
        sp.call(['cp',CONTCAR,'POSCAR'])
        if runType == 'relaxed':
            sp.call('cp INCAR.isif2 INCAR'.split())
        else:
            sp.call('cp INCAR.static INCAR'.split())
        # create submission script, place copy in appropriate subdirectory
        genSubScript(cname,tList,runLength,NCORES)
        # sp.call(('cp submission_script '+cname+"_%.5f/"%t).split())
        cell = Cell().loadFromPOSCAR('POSCAR')
        # strain tensors from Cerny 2004
        if 'cprime' in cname:
            strainCell(cell, [t,-t,0,0,0,0])
        elif 'c44' in cname:
            strainCell(cell, [0,0,0,t,0,0])
        else: # 'c66' in cname
            strainCell(cell, [0,0,0,0,0,t])
        cell.setSiteMobilities(True,True,True) # allow all directions to relax
        cell.sendToPOSCAR()
        # copy files to subdirectory, move subdirectory to $WORK
        sp.call(['mkdir',cname+'_%.5f'%t])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+[cname+'_submit']+\
                [cname+'_%.5f/'%t])
        sp.call('cp -r '+cname+'_%.5f '%t + WORK,shell=True)
        sp.call('chmod u+x %s_submit'%cname,shell=True)
    # run submission script
    sp.call(['sbatch','%s_submit'%cname])

#===========================================================================
# MAIN PROGRAM
#===========================================================================
runType = 'unrelaxed' # relaxed or unrelaxed calculations
runLength = 5 # Maximum run time in minutes
PERC = 1.0 # Maximum percent strain
NT = 5 # Number of strains to run per parameter

tmax = PERC/100
tList = []
# create list of strain amounts to run
for x in range(1, NT/2 + 1):
    tList += [tmax/x,-tmax/x]
# if odd # of strains, add one extra to account for integer division
if (NT%2==1):
    tList.append(tmax/(NT/2+1))
tList.sort()
print tList

eList = []
for dir in os.listdir(HOME):
    dirInfo = dir.split('_')
    eList += dirInfo[0] # main strain amount

for e in eList:
    # run Cprime calculations
    getCxx('%.5f_cprime'%e, e, tList, runType, runLength, NCORES)
    # run C44 calculations
    getCxx('%.5f_c44'%e, e,tList, runType, runLength, NCORES)
    # run C66 calculations
    getCxx('%.5f_c66'%e, e, tList, runType, runLength, NCORES)