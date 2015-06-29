# Jonas Kaufman jlkaufman@gmail.com
# June 22 2014

import subprocess as sp
from Cell import *
import math
import numpy as np

# number of cores for each calculation
NCORES = 32
#run length in minutes
RUNLENGTH = 120
# Stampede allocation number
ALLOCATION = 'TG-DMR140093'
# home and work directories
HOME = '/home1/03022/bassman/Jonas/'
#WORK = '/work/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'
# email address for Slurm notifications
EMAIL = 'jlkaufman@hmc.edu'

# functions to add layers
def addA(cell,z):
    cell.addSite(0,Site([0,0,z]))
    cell.addSite(0,Site([0.5,0.5,z]))
def addB(cell,z):
    cell.addSite(0,Site([0,(2.0/6.0),z]))
    cell.addSite(0,Site([0.5,(5.0/6.0),z]))
def addC(cell,z):
    cell.addSite(0,Site([0,(4.0/6.0),z]))
    cell.addSite(0,Site([0.5,(1.0/6.0),z]))

# add ESF and ISF functions
def makeISF(nshifts,nlayers,ngap):
    halfz = .5*(nlayers-1)/(nlayers+ngap) # z position right between the surface layers
    toptwo = 1.5/(nlayers+ngap)
    bottomtwo = (nlayers-2.5)/(nlayers+ngap)
    num = 0
    displacements = np.linspace(0.0, 1.0, nshifts+1)
    for disp in displacements:    
        cell = Cell().loadFromPOSCAR()
        cell.setSiteMobilities(False,False,True)
        for element in cell.sites:
            for i in range(cell.numberOfElements()):
                for j in range(cell.numberOfAtomsOfElement(i)):
                        x,y,z = cell.sites[i][j].position
                        if z > halfz:
                            cell.sites[i][j].move([x, y - disp*(1.0/6.0), z])
                        if z < toptwo or z > bottomtwo:
                            cell.sites[i][j].zfree = False
        cell.sendToPOSCAR('POSCAR_%02d'%num)
        num += 1

def makeESF(nshifts,nlayers,ngap):
    toptwo = 1.5/(nlayers+ngap)
    bottomtwo = (nlayers-2.5)/(nlayers+ngap)
    halfz = .5*(nlayers-3)/(nlayers+ngap) # z position below the ISF
    num = nshifts+1
    displacements = np.linspace(0.0, 1.0, num)
    for disp in displacements[1:]:    
        cell = Cell().loadFromPOSCAR('POSCAR_%02d'%(nshifts))
        cell.setSiteMobilities(False,False,True)
        for element in cell.sites:
            for i in range(cell.numberOfElements()):
                for j in range(cell.numberOfAtomsOfElement(i)):
                        x,y,z = cell.sites[i][j].position
                        if z < halfz:
                            cell.sites[i][j].move([x, y + disp*(1.0/6.0), z])
                        if z < toptwo or z > bottomtwo:
                            cell.sites[i][j].zfree = False
        cell.sendToPOSCAR('POSCAR_%02d'%num)
        num += 1

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
    '#SBATCH --mail-type=all\n' +               # send all emails
    '#SBATCH -A ' + ALLOCATION + '\n' +         # specifies project
    'module load vasp\n')                       # load vasp module
    for i in range(len(aList)):
        # change to work directory, run vasp
        string += 'cd '+WORK+'%2d_%s\n'%(aList[i],cname)
        string += "ibrun -o %d -n %d vasp_std > vasp_output.out &\n"%(NCORES*i,NCORES)
    string += 'wait\ncd '+HOME+'\nmkdir %s_results\n'%(cname)
    for i in range(len(aList)):
        # move directories to results directory
        string += 'cd '+WORK+'\nmv %2d_%s '%(aList[i],cname)
        string += HOME+'%s_results/\n'%cname
    f = open(cname + '_submit','w')
    f.write(string)
    f.close()  

def runJobs(cname, aList,runLength,NCORES):
    """
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK/Jonas directory and runs submission scripts
    """
    for a in aList:
        # create submission script
        genSubScript(cname,aList,runLength,NCORES)
        # copy files to subdirectory, move subdirectory to WORK
        sp.call(('cp POSCAR_%02d POSCAR'%a).split())
        sp.call(['mkdir','%02d_%s'%(a,cname)])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+[cname+'_submit']+\
                ['%02d_%s'%(a,cname)])
        sp.call('cp -r %02d_%s '%(a,cname)+WORK,shell=True)
        sp.call('chmod u+x %s_submit'%cname,shell=True)
    # run submission script
    sp.call(['sbatch','%s_submit'%cname])  
#===========================================================================
# MAIN PROGRAM
#===========================================================================
"""
element = 'Cu'
a = 3.63548
nshifts = 5
jobName = 'Cu_GSF'
"""
ngap = 3

element = raw_input('Element: ')
a = float(raw_input('Lattice parameter: '))
# ngap = int(raw_input('Vacuum gap size (number of atoms): '))
nshifts = int(raw_input('Number of points: '))
jobName = raw_input('Job name: ')

sequence = 'ABCABCABCABC'
layers = list(sequence)
nlayers = len(layers)
ax = math.sqrt(2)/2
ay = math.sqrt(6)/2
az = math.sqrt(3)/3
bx,by,bz = ax,ay,az*(nlayers+ngap)

# make perfect crystal to start from
cell = Cell() # Cartesian coordinates by default
cell.setHeader(element+' fcc')
cell.setElements([element]) # fix this?
cell.setA0(a)
cell.setCoordinateSystem('Direct')
cell.setLatticeVectors([[bx,0,0],[0,by,0],[0,0,bz]])
zpoint = 0.0
while layers:
    next = layers.pop(0)
    if next == 'A':
        addA(cell,zpoint)
    elif next == 'B':
        addB(cell,zpoint)
    else: # next == 'C'
        addC(cell,zpoint)
    zpoint += (1.0/(nlayers+ngap))
cell.setSiteMobilities(False,False,True)
cell.sendToPOSCAR()

makeISF(nshifts,nlayers,ngap) # add ISF
makeESF(nshifts,nlayers,ngap) # add ESF

aList = range(2*nshifts+1)
runJobs(jobName,aList,RUNLENGTH,NCORES) # run jobs