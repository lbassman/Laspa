#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  June 22 2014
#==============================================================================
"""
Add POTCAR, KPOINTS, and INCAR files to the working directory
Make script executable using 'chmod +x _____.py' to call as bash script
"""
import subprocess as sp
from Cell import *
import math
import numpy as np
# home and work directories (SET THESE TO YOUR OWN)
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'
# email address for Slurm notifications (SET TO YOUR OWN)
EMAIL = 'jlkaufman@hmc.edu'
# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

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
def makeISF(nShifts,nLayers,ny,nGap):
    halfz = .5*(nLayers-1)/(nLayers+nGap) # z position right between the surface layers
    topTwo = 1.5/(nLayers+nGap)
    bottomTwo = (nLayers-2.5)/(nLayers+nGap)
    displacements = np.linspace(0.0, 1.0, nShifts+1)
    for disp in displacements:    
        cell = Cell().loadFromPOSCAR()
        cell.setSiteMobilities(False,False,True)
        for element in cell.sites:
            for i in range(cell.numberOfElements()):
                for j in range(cell.numberOfAtomsOfElement(i)):
                        x,y,z = cell.sites[i][j].position
                        if z > halfz:
                            #cell.sites[i][j].move([x-disp*0.5, y - disp*(1.0/6.0)/ny, z])
                            cell.sites[i][j].move([x, y + disp*(1.0/3.0)/ny, z])
                        if z < topTwo or z > bottomTwo:
                            cell.sites[i][j].zfree = False
        cell.sendToPOSCAR(('POSCAR_%.5f'%disp).replace('.',''))

def makeESF(nShifts,nLayers,ny,nGap):
    topTwo = 1.5/(nLayers+nGap)
    bottomTwo = (nLayers-2.5)/(nLayers+nGap)
    halfz = .5*(nLayers-3)/(nLayers+nGap) # z position below the ISF
    displacements = np.linspace(0.0, 1.0, nShifts+1)
    for disp in displacements[1:]:    
        cell = Cell().loadFromPOSCAR(('POSCAR_%.5f'%(1.0)).replace('.',''))
        cell.setSiteMobilities(False,False,True)
        for element in cell.sites:
            for i in range(cell.numberOfElements()):
                for j in range(cell.numberOfAtomsOfElement(i)):
                        x,y,z = cell.sites[i][j].position
                        if z < halfz:
                            #cell.sites[i][j].move([x+disp*0.5, y+disp*(1.0/6.0)/ny, z])
                            cell.sites[i][j].move([x, y+disp*(1.0/3.0)/ny, z])
                        if z < topTwo or z > bottomTwo:
                            cell.sites[i][j].zfree = False
        num = disp + 1.0
        cell.sendToPOSCAR(('POSCAR_%.5f'%num).replace('.',''))

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

def runJobs(jName, aList,runLength,NCORES):
    """
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK/Jonas directory and runs submission scripts
    """
    dirList = []
    for a in aList:
        # copy files to subdirectory, move subdirectory to WORK
        dirName = '%s_%.5f'%(jName,a)
        dirList += [dirName]
        sp.call((('cp POSCAR_%.5f POSCAR'%a).replace('.','')).split())
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
elements = raw_input('Element: ')
if not elements: elements = ['Cu']
else: elements = [elements]
print ''.join(elements),'\n'
a = raw_input('Lattice constant (angstroms): ')
if not a: a= 3.6355
else: a = float(a)
print a,'\n'
nGap = raw_input('Vacuum gap width (number of atoms): ')
if not nGap: nGap = 4
else: nGap = int(nGap)
print nGap,'\n'
# periodicity in x
ny = raw_input('Periodicity in y: ')
if not ny: ny = 1
else: ny = int(ny)
print ny,'\n'
nShifts = raw_input('Number shifts per stacking fault: ')
if not nShifts: nShifts = 10
else: nShifts = int(nShifts)
print nShifts,'\n'
pathway = raw_input('Fault pathway (ISF only or ESF): ')
if 'e' in pathway or 'E' in pathway:
    pathway = 'ESF'
else:
    pathway = 'ISF'
print pathway,'\n'
runLength = raw_input('Maximum run time (minutes): ')
if not runLength: runLength = 300
else: runLength = int(runLength)
print runLength,'\n'
jobName = raw_input('Job name: ')
if not jobName: jobName = 'GSF'
print jobName,'\n'
NCORES = raw_input('Number of cores per simulation: ')
if not NCORES: NCORES = 16
else: NCORES = int(NCORES)
print NCORES,'\n'
resultsDir = raw_input('Put results in home or work: ')
if 'w' in resultsDir or 'W' in resultsDir: HOME = WORK
print HOME,'\n'

sequence = 'ABCABCABCABC'
layers = list(sequence)
nLayers = len(layers)
ax = math.sqrt(2)/2
ay = math.sqrt(6)/2
az = math.sqrt(3)/3
bx,by,bz = ax,ny*ay,az*(nLayers+nGap)

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
    zpoint += (1.0/(nLayers+nGap))
cell.setSiteMobilities(False,False,True)
cell.sendToPOSCAR()

totalDisp = 1.0
points = nShifts + 1
makeISF(nShifts,nLayers,ny,nGap) # add ISF
if pathway == 'ESF':
    makeESF(nShifts,nLayers,ny,nGap) # add ESF
    totalDisp = 2.0
    points = 2*nShifts+1
aList = np.linspace(0.0, totalDisp,points)
runJobs(jobName,aList,runLength,NCORES) # run jobs