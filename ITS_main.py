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

#SCRIPT = 'IYS_main.py'

# number of cores for each calculation
NCORES = 16 ## should this be 32?
runLength = 5 ## to be set by user later on

# home and work directories
HOME = '/home1/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'

# email address for Slurm notifications
EMAIL = 'jlkaufman@hmc.edu'

# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

# default increment and max strain, decimal
INC = 0.01
EMAX = 0.02

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

def writeState(eN, inc, emax, message=''):
    # update state file and close
    #decimals, not percents
    state = open('STATE','w')
    string = 'eN '+str(eN)+'\ninc '+str(inc)+'\nemax '+str(emax)+message
    state.write(string)
    state.close()
    return

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
    #'#SBATCH -v\n' +                            # verbose 
    '#SBATCH -J ' + cname +  '\n' +             # specify job name
    '#SBATCH -o ' + cname + '%j\n' +            # write output to this file
    '#SBATCH -n %d\n'%(NCORES*len(tList)) +     # request cores
    '#SBATCH -p normal\n' +                     # send to normal queue
    '#SBATCH -t %.2d:%d:00\n'%(hrs,mins) +      # set maximum wall (clock) time
    '#SBATCH --mail-user=' + EMAIL +'\n' +      # set email
    #'#SBATCH --mail-type=begin\n' +             # send email when job starts
    #'#SBATCH --mail-type=end\n' +               # send email when job ends
    '#SBATCH --mail-type=all\n' +               
    '#SBATCH -A ' + ALLOCATION + '\n' +         # specifies project
    'module load vasp\n')                       # load vasp module
    for i in range(len(tList)):
        # change to work directory, run vasp
        string += 'cd '+WORK+'%.5f_%s\n'%(tList[i],cname)
        string += "ibrun -o %d -n %d vasp_std > vasp_output.out &\n"%(NCORES*i,NCORES)
    string += 'wait\ncd '+HOME+'\nmkdir %.5f_%s_results\n'%(tList[i],cname)
    for i in range(len(tList)):
        # move directories to results directory
        string += 'cd '+WORK+'\nmv %.5f_%s '%(tList[i],cname)
        string += HOME+'%.5f_%s_results\n'%(tList[i],cname)
    f = open(cname + '_submit','w')
    f.write(string)
    f.close()

def relaxCell(jobName, eN,runLength,NCORES):
    """
    generate subdirectory, copies files, and submits job to
    relax a strained cell
    """
    tList = [eN]
    cname = jobName
    genSubScript(cname,tList,runLength,NCORES)
    # copy files to subdirectory, move subdirectory to WORK
    for t in tList:
        # create submission script
        genSubScript(cname,tList,runLength,NCORES)
        # copy files to subdirectory, move subdirectory to WORK
        sp.call(['mkdir','%.5f_%s'%(t,cname)])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+[cname+'_submit']+\
                ['%.5f_%s'%(t,cname)])
        sp.call('cp -r %.5f_%s '%(t,cname)+WORK,shell=True)
        sp.call('chmod u+x %s_submit'%cname,shell=True)
    sp.call(['sbatch','%s_submit'%cname])

#===========================================================================
# MAIN PROGRAM
#===========================================================================
# user input, to be added later
"""
if no STATE FILE, ask for user input, otherwise get information from previous run

NCORES = float(raw_input('Number of cores per simulation: '))
emax = float(raw_input('Maximum strain (%): '))/100.0     # e.g. 15%
inc = float(raw_input('Strain step size (%): '))/100.0    # e.g. 1%
numSteps = int(emax/inc)
valid = 0
while not valid:
    try:
        runLength = int(raw_input("Maximum run time in minutes: "))
        valid = 1
    except:
        print "Run time must be a whole number"

# start with relaxed cell
sp.call('cp CONTCAR.E0 POSCAR'.split())
sp.call('cp CONTCAR.E0 CONTCAR'.split())
"""
#######################################################################
cname = 'main'
# read in last state
state = readState()
eN = state[0]
inc = state[1]
emax = state[2]
aN = eN + 1
if abs(eN) < (abs(emax)+inc/2.0):
    print '\nRunning next strain step.\n'
    # strain cell by increment relative to starting position
    # a(n+1) = f(n)*a(n), f(n) = 1 + inc*a(0)/a(n)
    a0 = aN
    f = 1+inc/a0
    if eN == 0.0:
        CONTCAR = 'CONTCAR.E0'
    else:
        CONTCAR = HOME+'%.5f_%s_results/'%(eN,cname)+'%.5f_%s/CONTCAR'%(eN,cname)
    cell = Cell().loadFromPOSCAR(CONTCAR) # use previous CONTCAR as POSCAR
    strainCell(cell,[f-1,0,0,0,0,0])
    cell.setSiteMobilities(False,True,True)
    cell.sendToPOSCAR()     
    aN = a0*f
    eN = aN -1 # this is the next strain step
    
    # write out final state
    writeState(eN, inc, emax)
    # run VASP relaxation
    relaxCell(cname,eN,runLength,NCORES)
else:
    writeState(eN,inc,emax,'\nReached maximum strain.\n')
    sys.exit('\nReached maximum strain.\n')