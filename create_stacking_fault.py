#!/usr/bin/env python

# Jonas Kaufman jlkaufman@gmail.com
# June 22 2014

import subprocess as sp
from Cell import *
import math
import numpy as np

# Stampede allocation number
ALLOCATION = 'TG-DMR140093'
# home and work directories
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'
# email address for Slurm notifications
EMAIL = 'jlkaufman@hmc.edu'

# functions to add layers
def addA(cell,ny,z):
    y = 0.0
    while y < .99: # .99 instead of 1.0 to account for floating point error
        cell.addSite(0,Site([0.0,y,z]))
        y += 1.0/(2*ny)
        cell.addSite(0,Site([0.5,y,z]))
        y += 1.0/(2*ny)
def addB(cell,ny,z):
    y = 2.0/(6.0*ny)
    while y < .99:
        cell.addSite(0,Site([0.0,y,z]))
        y += 1.0/(2*ny)
        cell.addSite(0,Site([0.5,y,z]))
        y += 1.0/(2*ny)
def addC(cell,ny,z):
    y = 1.0/(6.0*ny)
    while y < .99:
        cell.addSite(0,Site([0.5,y,z]))
        y += 1.0/(2*ny)
        cell.addSite(0,Site([0.0,y,z]))
        y += 1.0/(2*ny)

# add ESF and ISF functions
def makeISF(nshifts,nlayers,ny,ngap):
    halfz = .5*(nlayers-1)/(nlayers+ngap) # z position right between the surface layers
    toptwo = 1.5/(nlayers+ngap)
    bottomtwo = (nlayers-2.5)/(nlayers+ngap)
    displacements = np.linspace(0.0, 1.0, nshifts+1)
    for disp in displacements:    
        cell = Cell().loadFromPOSCAR()
        cell.setSiteMobilities(False,False,True)
        for element in cell.sites:
            for i in range(cell.numberOfElements()):
                for j in range(cell.numberOfAtomsOfElement(i)):
                        x,y,z = cell.sites[i][j].position
                        if z > halfz:
                            cell.sites[i][j].move([x-disp*0.5, y - disp*(1.0/6.0)/ny, z])
                        if z < toptwo or z > bottomtwo:
                            cell.sites[i][j].zfree = False
        cell.sendToPOSCAR(('POSCAR_%.5f'%disp).replace('.',''))

def makeESF(nshifts,nlayers,ny,ngap):
    toptwo = 1.5/(nlayers+ngap)
    bottomtwo = (nlayers-2.5)/(nlayers+ngap)
    halfz = .5*(nlayers-3)/(nlayers+ngap) # z position below the ISF
    displacements = np.linspace(0.0, 1.0, nshifts+1)
    for disp in displacements[1:]:    
        cell = Cell().loadFromPOSCAR(('POSCAR_%.5f'%(1.0)).replace('.',''))
        cell.setSiteMobilities(False,False,True)
        for element in cell.sites:
            for i in range(cell.numberOfElements()):
                for j in range(cell.numberOfAtomsOfElement(i)):
                        x,y,z = cell.sites[i][j].position
                        if z < halfz:
                            cell.sites[i][j].move([x+disp*0.5, y+disp*(1.0/6.0)/ny, z])
                        if z < toptwo or z > bottomtwo:
                            cell.sites[i][j].zfree = False
        num = disp + 1.0
        cell.sendToPOSCAR(('POSCAR_%.5f'%num).replace('.',''))

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
    '#SBATCH -n %d\n'%(NCORES*len(aList)) +     # request 64 cores
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

def runJobs(jName, aList,runLength,NCORES):
    """
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK/Jonas directory and runs submission scripts
    """
    for a in aList:
        # create submission script
        genSubScript(jName,aList,runLength,NCORES)
        # copy files to subdirectory, move subdirectory to WORK
        sp.call((('cp POSCAR_%.5f POSCAR'%a).replace('.','')).split())
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
# USER INPUTS
# which pure element is in system
elements = raw_input('Element: ')
if not elements: elements = ['Cu']
else: elements = [elements]
print ''.join(elements)
# lattice constant
a = raw_input('Lattice constant (angstroms): ')
if not a: a= 3.6355
else: a = float(a)
print a
# width of vaccuum gap
ngap = raw_input('Vacuum gap width (number of atoms): ')
if not ngap: ngap = 4
else: ngap = int(ngap)
print ngap
# periodicity in x
###########################
# periodicity in y
ny = raw_input('Periodicity in y: ')
if not ny: ny = 1
else: ny = int(ny)
print ny
# number of shifts per SF
nshifts = raw_input('Number shifts per stacking fault: ')
if not nshifts: nshifts = 10
else: nshifts = int(nshifts)
print nshifts
# pathway (ISF or full)
###########################
# run time
runLength = raw_input('Maximum run time (minutes): ')
if not runLength: runLength = 300
else: runLength = int(runLength)
print runLength
# job name
jobName = raw_input('Job name: ')
if not jobName: jobName = 'GSF'
print jobName
# number of cores
NCORES = raw_input('Number of cores per simulation: ')
if not NCORES: NCORES = 16
else: NCORES = int(NCORES)
print NCORES

sequence = 'ABCABCABCABC'
layers = list(sequence)
nlayers = len(layers)
ax = math.sqrt(2)/2
ay = math.sqrt(6)/2
az = math.sqrt(3)/3
bx,by,bz = ax,ny*ay,az*(nlayers+ngap)

# make perfect crystal to start from
cell = Cell() # Cartesian coordinates by default
cell.setHeader(jobName)
cell.setElements(elements)
cell.setA0(a)
cell.setCoordinateSystem('Direct') # change to direct
cell.setLatticeVectors([[bx,0,0],[0,by,0],[0,0,bz]])
zpoint = 0.0
while layers:
    next = layers.pop(0)
    if next == 'A':
        addA(cell,ny,zpoint)
    elif next == 'B':
        addB(cell,ny,zpoint)
    else: # next == 'C'
        addC(cell,ny,zpoint)
    zpoint += (1.0/(nlayers+ngap))
cell.setSiteMobilities(False,False,True)
cell.sendToPOSCAR()

makeISF(nshifts,nlayers,ny,ngap) # add ISF
makeESF(nshifts,nlayers,ny,ngap) # add ESF

aList = np.linspace(0.0, 2.0, 2*nshifts+1)
runJobs(jobName,aList,runLength,NCORES) # run jobs