# Jonas Kaufman jlkaufman@hmc.edu
# February 25, 2015
# Script to run VASP calculations necessary to find lattice parameters
# of pure element in fcc, bcc (or hcp) crystal structure - on Stampede
# Based on a scipt by Josh Sanz

# this script assumes these lattice vectors:
#  1.00000         0.00000000000000      0.00000
# -0.50000         0.86602540378444      0.00000
#  0.00000         0.00000000000000      c/a
# and that a is the scaling factor

"""
Add the POSCAR, POTCAR, KPOINTS, and INCAR files to the working directory.

Run this script to generate the necessary POSCAR files and submit all jobs
to the queue on Stampede.
"""

from Cell import *
import subprocess as sp
import numpy as np

# number of cores for each calculation
NCORES = 16

# home and work directories
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'

# email address for Slurm notifications
EMAIL = 'jlkaufman@hmc.edu'

# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

def genSubScript(jName,dirList,runLength,NCORES):
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
    '#SBATCH -n %d\n'%(NCORES*len(dirList)) +   # request cores
    '#SBATCH -p normal\n' +                     # send to normal queue
    '#SBATCH -t %02d:%02d:00\n'%(hrs,mins) +    # set maximum wall (clock) time
    '#SBATCH --mail-user=' + EMAIL +'\n' +      # set email
    '#SBATCH --mail-type=all\n' +               # send all emails
    '#SBATCH -A ' + ALLOCATION + '\n' +         # specifies project
    'module load vasp\n')                       # load vasp module
    for i in range(len(dirList)):
        # change to work directory, run vasp
        string += 'cd '+WORK+'%s\n'%dirList[i]
        string += "ibrun -o %d -n %d vasp_std > vasp_output.out &\n"%(NCORES*i,NCORES)
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
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK/Jonas directory and runs submission scripts
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
    # create submission script
    genSubScript(jName,dirList,runLength,NCORES)
    sp.call('chmod u+x %s_submit'%jName,shell=True)
    # run submission script
    sp.call(['sbatch','%s_submit'%jName]) 
def getLatHex(jName, aList, caList, runLength,NCORES):
    """
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK/Jonas directory and runs submission scripts
    """
    dirList = []
    for a in aList:
        for ca in caList:
            cell = Cell().loadFromPOSCAR()
            cell.setA0(float(a))
            (x,y,z) = cell.latticeVectors
            # this script assumes these lattice vectors:
            #  1.00000         0.00000000000000      0.00000
            # -0.50000         0.86602540378444      0.00000
            #  0.00000         0.00000000000000      c/a
            # and that a is the scaling factor
            cell.setLatticeVectors([x,y,[0,0,ca]])
            cell.sendToPOSCAR()    
            # copy files to subdirectory, move subdirectory to WORK
            dirName = '%s_%.5f_%.5f'%(jName,a,ca)
            dirList += [dirName]
            sp.call(['mkdir',dirName])
            #sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+[jName+'_submit']+\
            #        [dirName])
            sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+\
                    [dirName])
            sp.call(('cp -r %s '%dirName)+WORK,shell=True)
    # create submission script (not sure why this was being done for each
    genSubScript(jName,dirList,runLength,NCORES)  
    sp.call('chmod u+x %s_submit'%jName,shell=True)
    sp.call(['sbatch','%s_submit'%jName])   

hcp = False
# add default values
runTime = int(raw_input('Maximum run time in minutes: '))
structure = raw_input('Crystal structure (fcc, bcc or hcp): ')
if 'h' in structure or 'H' in structure:
    hcp = True
    ca_min = float(raw_input('Minimum c/a ratio: '))
    ca_max = float(raw_input('Maximum c/a ratio: '))
    nca_points = float(raw_input('Number of c/a values: '))
a_min = float(raw_input('Minimum lattice parameter a in angstroms: '))
a_max = float(raw_input('Maximum lattice parameter a in angstroms: '))
na_points = int(raw_input('Number of a values: '))
jName = raw_input('Job name: ')
resultsDir = raw_input('Results in home or work: ')
if 'w' in resultsDir or 'W' in resultsDir: HOME = WORK
# make list of lattice parameters
aList = np.linspace(a_min, a_max, na_points)
print aList
if hcp:
    caList = np.linspace(ca_min, ca_max, nca_points) 
    print caList
    getLatHex(jName,aList,caList,runTime,NCORES)
else:
    getLatCubic(jName,aList,runTime,NCORES)
