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

def genSubScript(cname,aList,runLength,NCORES):
    """
    create a submission script for Stampede's SLURM queueing system
    this version sends an email whenever one of the queued jobs starts
    """

    # uses integer division
    hrs = runLength/60
    mins = runLength%60
    string = ('#!/bin/bash\n' +
    '#SBATCH -J ' + cname +  '\n' +             # specify job name
    '#SBATCH -o ' + cname + '%j\n' +            # write output to this file
    '#SBATCH -n %d\n'%(NCORES*len(aList)) +     # request 64 cores
    '#SBATCH -p normal\n' +                     # send to normal queue
    '#SBATCH -t %.2d:%d:00\n'%(hrs,mins) +      # set maximum wall (clock) time
    '#SBATCH --mail-user=' + EMAIL +'\n' +      # set email
    '#SBATCH --mail-type=begin\n' +             # send email when job starts
    '#SBATCH --mail-type=end\n' +               # send email when job ends
    '#SBATCH -A ' + ALLOCATION + '\n' +         # specifies project
    'module load vasp\n')                       # load vasp module
    for i in range(len(aList)):
        # change to work directory, run vasp
        string += 'cd '+WORK+'%.5f_%s\n'%(aList[i],cname)
        string += "ibrun -o %d -n %d vasp_std > vasp_output.out &\n"%(NCORES*i,NCORES)
    string += 'wait\ncd '+HOME+'\nmkdir %s_results\n'%(cname)
    for i in range(len(aList)):
        # move directories to results directory
        string += 'cd '+WORK+'\nmv %.5f_%s '%(aList[i],cname)
        string += HOME+'%s_results/\n'%cname
    f = open(cname + '_submit','w')
    f.write(string)
    f.close()

def getLat(cname, aList,runLength,NCORES):
    """
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK/Jonas directory and runs submission scripts
    """
    for a in aList:
        cell = Cell().loadFromPOSCAR('POSCAR')
        cell.setA0(float(a))
        cell.sendToPOSCAR()
        # create submission script
        genSubScript(cname,aList,runLength,NCORES)
        # copy files to subdirectory, move subdirectory to WORK
        sp.call(['mkdir','%.5f_%s'%(a,cname)])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+[cname+'_submit']+\
                ['%.5f_%s'%(a,cname)])
        sp.call('cp -r %.5f_%s '%(a,cname)+WORK,shell=True)
        sp.call('chmod u+x %s_submit'%cname,shell=True)
    # run submission script
    sp.call(['sbatch','%s_submit'%cname])    


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
cname = raw_input('Job name: ')

# make list of lattice parameters
aList = np.linspace(a_min, a_max, n_points)
print aList
# run simulations
getLat(cname,aList,runTime,NCORES)