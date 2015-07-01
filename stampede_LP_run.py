# Jonas Kaufman jlkaufman@hmc.edu
# February 25, 2015
# Script to run VASP calculations necessary to find lattice parameters
# of pure element in fcc, bcc (or hcp) crystal structure - on Stampede
# Based on a scipt by Josh Sanz

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
HOME = '/home1/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'

# email address for Slurm notifications
EMAIL = 'jlkaufman@hmc.edu'

# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

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

def getLat(jName, aList,runLength,NCORES):
    """
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK/Jonas directory and runs submission scripts
    """
    for a in aList:
        cell = Cell().loadFromPOSCAR('POSCAR')
        cell.setA0(float(a))
        cell.sendToPOSCAR()
        # create submission script
        genSubScript(jName,aList,runLength,NCORES)
        # copy files to subdirectory, move subdirectory to WORK
        sp.call(['mkdir','%s_%.5f'%(jName,a)])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+[jName+'_submit']+\
                ['%s_%.5f'%(jName,a)])
        sp.call('cp -r %s_%.5f '%(jName,a)+WORK,shell=True)
        sp.call('chmod u+x %s_submit'%jName,shell=True)
    # run submission script
    sp.call(['sbatch','%s_submit'%jName])    


valid = 0
while not valid:
    try:
        runTime = int(raw_input("Maximum run time in minutes: "))
        valid = 1
    except:
        print "Run time must be a whole number"

a_min = float(raw_input("Minimum lattice parameter in angstroms: "))
a_max = float(raw_input("Maximum lattice parameter in angstroms: "))
n_points = int(raw_input("Number of points: "))
jName = raw_input('Job name: ')

# make list of lattice parameters
aList = np.linspace(a_min, a_max, n_points)
print aList
# run simulations
getLat(jName,aList,runTime,NCORES)