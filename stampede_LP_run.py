#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  February 25, 2015
#  Script to run VASP calculations necessary to find lattice parameters
#  of a fcc, bcc or hcp crystal structure - on Stampede
#==============================================================================
# This script assumes these lattice vectors for hcp:
#  1.00000         0.00000000000000      0.00000
# -0.50000         0.86602540378444      0.00000
#  0.00000         0.00000000000000      c/a
# where a is the overall scaling factor
"""
Add POSCAR, POTCAR, KPOINTS, and INCAR files to the working directory
Make script executable using 'chmod +x _____.py' to call as bash script
"""
from Cell import *
import subprocess as sp
import numpy as np
# home and work directories (SET THESE TO YOUR OWN)
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'
# email address for Slurm notifications (SET TO YOUR OWN)
EMAIL = 'jlkaufman@hmc.edu'
# Stampede allocation number
ALLOCATION = 'TG-DMR140093'
# number of cores for each calculation (16 should be fine)
NCORES = 16

def genSubScript(jName,dirList,runLength,NCORES):
    """ creates a submission script for Stampede's SLURM queueing system """
    # uses integer division
    hrs = runLength/60
    mins = runLength%60
    string = ('#!/bin/bash\n' +
    '#SBATCH -J ' + jName +  '\n' +           # specify job name
    '#SBATCH -o ' + jName + '%j\n' +          # write output to this file
    '#SBATCH -n %d\n'%(NCORES*len(dirList)) + # request cores
    '#SBATCH -p normal\n' +                   # send to normal queue
    '#SBATCH -t %02d:%02d:00\n'%(hrs,mins) +  # set maximum wall time
    '#SBATCH --mail-user=' + EMAIL +'\n' +    # set email
    '#SBATCH --mail-type=all\n' +             # send all emails
    '#SBATCH -A ' + ALLOCATION + '\n' +       # specify project
    'module load vasp\n')                     # load vasp module
    for i in range(len(dirList)):
        # change to work directory, run vasp
        string += 'cd '+WORK+'%s\n'%dirList[i]
        string += "ibrun -o %d -n %d vasp_std > vasp_output.out &\n"%(NCORES*i,NCORES)
    # wait for all jobs to finish, move to results directory
    string += 'wait\ncd '+HOME+'\nmkdir %s_results\n'%(jName)
    for i in range(len(dirList)):
        # move directories to results directory
        string += 'cd '+WORK+'\nmv %s '%dirList[i]
        string += HOME+'%s_results/\n'%jName
    f = open(jName + '_submit','w')
    f.write(string)
    f.close()

def getLatCubic(jName, aList,runLength,NCORES):
    """
    creates the necessary POSCARs for a cubic structure
    generates a subdirectory for each vasp run, each with the necessary files,
    moves subdirectories to work directory and runs submission script
    """
    dirList = []
    for a in aList:
        cell = Cell().loadFromPOSCAR()
        cell.setA0(float(a))
        cell.sendToPOSCAR()
        # copy files to subdirectory, move subdirectory to WORK
        dirName = '%s_%.5f'%(jName,a)
        dirList += [dirName]
        sp.call(['mkdir',dirName])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+\
                [dirName])
        sp.call(('cp -r %s '%dirName)+WORK,shell=True)
    # create submission script and run
    genSubScript(jName,dirList,runLength,NCORES)
    sp.call('chmod u+x %s_submit'%jName,shell=True)
    sp.call(['sbatch','%s_submit'%jName]) 

def getLatHex(jName, aList, caList, runLength,NCORES):
    """
    creates the necessary POSCARs for a hexagonal structure
    generates a subdirectory for each vasp run, each with the necessary files,
    moves subdirectories to work directory and runs submission script
    """
    dirList = []
    for a in aList:
        for ca in caList:
            cell = Cell().loadFromPOSCAR()
            cell.setA0(float(a))
            (x,y,z) = cell.latticeVectors
            # assumes these lattice vectors:
            #  1.00000         0.00000000000000      0.00000
            # -0.50000         0.86602540378444      0.00000
            #  0.00000         0.00000000000000      c/a
            # where a is the overall scaling factor
            cell.setLatticeVectors([x,y,[0,0,ca]])
            cell.sendToPOSCAR()    
            # copy files to subdirectory, move subdirectory to WORK
            dirName = '%s_%.5f_%.5f'%(jName,a,ca)
            dirList += [dirName]
            sp.call(['mkdir',dirName])
            sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+\
                    [dirName])
            sp.call(('cp -r %s '%dirName)+WORK,shell=True)
    # create submission script and run
    genSubScript(jName,dirList,runLength,NCORES)  
    sp.call('chmod u+x %s_submit'%jName,shell=True)
    sp.call(['sbatch','%s_submit'%jName])   

#==============================================================================
#  Main Program
#==============================================================================
# get user inputs (with defaults)
structure = raw_input('Crystal structure (fcc, bcc or hcp): ')
hcp = False
if 'h' in structure or 'H' in structure:
    hcp = True
    structure = 'hcp'
elif 'b' in structure or 'B' in structure:
    structure = 'bcc'
else:
    structure = 'fcc'
print structure,'\n'
if hcp:
    caMin = raw_input('Minimum c/a ratio: ')
    if not caMin: caMin = 1.5
    else: caMin = float(caMin)
    print caMin,'\n'
    caMax = raw_input('Maximum c/a ratio: ')
    if not caMax: caMax = 2.0
    else: caMax = float(caMax)
    print caMax,'\n'
    caPoints = raw_input('Number of c/a values: ')
    if not caPoints: caPoints = 7
    else: caPoints = int(caPoints)
    print caPoints,'\n'  
aMin = raw_input('Minimum lattice parameter a in angstroms: ')
if not aMin: aMin = 2.0
else: aMin = float(aMin)
print aMin,'\n'
aMax = raw_input('Maximum lattice parameter a in angstroms: ')
if not aMax: aMax = 4.0
else: aMax = float(aMax)
print aMax,'\n'
aPoints = raw_input('Number of a values: ')
if not aPoints: aPoints = 7
else: aPoints = int(aPoints)
print aPoints,'\n'
runLength = raw_input('Maximum run time (minutes): ')
if not runLength: runLength = 300
else: runLength = int(runLength)
print runLength,'\n'
jName = raw_input('Job name: ')
if not jName: jName = 'LP'
print jName,'\n'
resultsDir = raw_input('Put results in home or work: ')
if 'w' in resultsDir or 'W' in resultsDir: HOME = WORK
print HOME,'\n'

# run jobs
aList = np.linspace(aMin, aMax, aPoints).tolist()
print 'a values:'
print aList,'\n'
if hcp:
    caList = np.linspace(caMin, caMax, caPoints).tolist()
    print 'c/a values:'
    print caList,'\n'
    getLatHex(jName,aList,caList,runLength,NCORES)
else:
    getLatCubic(jName,aList,runLength,NCORES)
