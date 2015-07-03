#!/usr/bin/env python

# Jonas Kaufman jlkaufman@hmc.edu
# June 29, 2015
# Script to parse VASP output directories located in
# parent directory: jobName_results

# Run 'module load python' before using on Stampede
# Make script executable using 'chmod +x _____.py'

import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
from pylab import * #this includes numpy as np!
from scipy.optimize import leastsq
import os
import fnmatch


# home and work directories (add $SCRATCH?)
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'

#===========================================================================
# Parsing Functions
#===========================================================================
def getTime(file):
    """ parses an OUTCAR file for the run time """
    f = open(file,'r')
    while True:
        nextLine = f.readline()
        if 'Total CPU time used' in nextLine:
            time = float(nextLine.split()[5])
            break
    return time

def getEnergies(file): # this function returns E(sigma->0), not TOTEN
    """ parses an OUTCAR file and pulls out the energy of the system
    after each ionic step """
    energies = []
    f = open(file,'r')
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' in nextLine:
            f.readline()    # line of dashes
            f.readline()    # TOTEN line
            f.readline()    # blank line
            energyLine = f.readline().split()
            energies += [float(energyLine[6])]
    return energies
    
def getPressures(file):
    """ parses an OUTCAR file and pulls out the pressure of the system
    after each ionic step """
    f = open(file,'r')
    pressures = []
    stresses = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'external pressure' in nextLine:
            PLine = nextLine.split()
            pressures += [float(PLine[3])]
            stresses += [float(PLine[8])]
    return (pressures,stresses)

def getSizes(file):
    """ returns the volume and lattice parameters after each ionic step """
    f = open(file,'r')
    volumes = []
    lats = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'VOLUME and BASIS-vectors are now :' in nextLine:
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            volumeLine = f.readline().split()
            volumes += [float(volumeLine[4])]
            for i in range(6):
                f.readline()    # text
            aLine = f.readline().split()
            ax = float(aLine[0])
            ay = float(aLine[1])
            az = float(aLine[2])
            lats.append([ax,ay,az])
    return (volumes, lats)    
    
def parseResults(directory):
    """ parses the files in each subdirectory of the given results directory """
    runs = []
    dList = []
    ELists = []
    VLists = []
    aLists = []
    PLists = []
    sLists = []
    tList = []
    for dir in os.listdir(directory):
        if os.path.isdir(directory + dir):
            OUTCAR = False
            print 'Found directory %s\nParsing...'%dir
            for file in os.listdir(directory + dir):        
                if fnmatch.fnmatch(file,'OUTCAR'):
                    pathToFile = directory+dir+'/'+file
                    dList += [dir] # add directory name
                    ELists += [getEnergies(pathToFile)] # add energy list
                    sizeData = getSizes(pathToFile) 
                    VLists += [sizeData[0]] # add volume list
                    aLists += [sizeData[1]] # add lattice parameters lists 
                    pressureData = getPressures(pathToFile)
                    PLists += [pressureData[0]] # add pressure lists
                    sLists += [pressureData[1]] # add Pullay stress lists
                    tList += [getTime(pathToFile)] # add time
                    print 'OUTCAR read\n'
                    OUTCAR = True
            if not OUTCAR:
                print 'WARNING: no OUTCAR read\n'
    for i in range(len(dList)):
        runs.append([dList[i],ELists[i],VLists[i],aLists[i],PLists[i],sLists[i],tList[i]])
    return runs

#===========================================================================
# Display and Organization Functions
#===========================================================================
def displayRun(run):
    dirName = run[0]
    EList = run[1]
    VList = run[2]
    aList = run[3]
    PList = run[4]
    sList = run[5]
    time = run[6]
    lines = []
    lines += [dirName]
    headings = ('E0','Volume','ax','ay','az','Pressure','Pullay')
    lines += [(('%-8s\t'*len(headings))%headings)] # need to fix columns
    # come up with better formatting
    nSteps = len(EList)
    for i in range(nSteps):
        lats = aList[i]
        ax = lats[0]
        ay = lats[1]
        az = lats[2]
        data = (EList[i],VList[i],ax,ay,az,PList[i],sList[i])
        lines += [(('%.6f\t'*len(data))%data)]
    lines += ['%d ionic steps'%nSteps]    
    lines += ['%d seconds'%time]
    print '\n'.join(lines)+'\n'

def displayFinal(runList):
    fins = finalValues(runList)
    lines = []
    lines += ['Final values']
    headings = ('E0','Volume','ax','ay','az','Pressure','Pullay')
    lines += [(('%-8s\t'*len(headings))%headings)]
    for i in range(len(runList)):
        E = fins[1][i]
        V = fins[2][i]
        lats = fins[3][i]
        ax = lats[0]
        ay = lats[1]
        az = lats[2]
        P = fins[4][i]
        s = fins[5][i]
        data = (E,V,ax,ay,az,P,s)
        lines += [(('%.6f\t'*len(data))%data)]
        # run time? dir name?
    print '\n'.join(lines)+'\n'

def finalValues(runList):
    dList = []
    EFins = []
    VFins = []
    aFins = []
    PFins = []
    sFins = []
    tList = []
    for run in runList:
        dList += [run[0]]
        EFins += [run[1][-1]]
        VFins += [run[2][-1]]
        aFins += [run[3][-1]]
        PFins += [run[4][-1]]
        sFins += [run[5][-1]]
        tList += [run[6]]
    return (dList,EFins,VFins,aFins,PFins,sFins,tList)

#===========================================================================
# Birch Murnaghan Fitting
#===========================================================================
def Birch(parameters,vol):
    '''
    given a vector of parameters and volumes, return a vector of energies.
    equation From Wikipedia
    '''
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    term12 = ((V0/vol)**(2.0/3.0) - 1.0)
    term3 = (6.0 - 4.0*(V0/vol)**(2.0/3.0))
    E = E0 + (9.0*V0*B0/16.0)*((term12**3.0)*BP + (term12**2.0)*term3)
    return E
# and we define an objective function that will be minimized
def objective(pars,y,x):
    #we will minimize this function
    err =  y - Birch(pars,x)
    return err

def fitBirch(EList,VList):
    v = np.array(VList)
    e = np.array(EList)
    vfit = np.linspace(min(v),max(v),100)
    ### fit a parabola to the data
    # y = ax^2 + bx + c
    a,b,c = polyfit(v,e,2) #this is from pylab
    #now here are our initial guesses.
    v0 = -b/(2*a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4

    x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Birch function

    birchpars, ier = leastsq(objective, x0, args=(e,v)) #this is from scipy
 
    #now we make a figure summarizing the results
    plot(v,e,'ro')
    plot(vfit, Birch(birchpars,vfit), label='Birch fit')
    xlabel('Volume ($\AA^3$)')
    ylabel('Energy (eV)')
    legend(loc='best')

    #add some text to the figure in figure coordinates
    ax = gca()
    text(0.4,0.7,'Min volume = %1.5f $\AA^3$' % birchpars[3],
         transform = ax.transAxes)
    text(0.4,0.6,'Min lattice constant = %1.5f $\AA$' % (birchpars[3]**(1.0/3.0)),
         transform = ax.transAxes)
    text(0.4,0.5,'Bulk modulus = %1.3f eV/$\AA^3$ = %1.3f GPa' % (birchpars[1],
                                                                  birchpars[1]*160.21773)
         , transform = ax.transAxes)


    savefig('birch_fit.png')
  
    """
    print 'initial guesses  : ',x0
    print 'fitted parameters: ', birchpars
    """

#===========================================================================
# MAIN PROGRAM
#===========================================================================
direct = raw_input("Results directory in home or work?: ")
if 'h' in direct or 'H' in direct:
    RESULTS = HOME
else:
    RESULTS = WORK
 # where all the results folders are located, change to be input
jName = raw_input('Job name: ') # check if this exists
runList = parseResults(RESULTS+jName+'_results/')
runList.sort(key=lambda x: x[0]) # sort job list by directory name

static = True
for run in runList:
    EList = run[1]
    if len(EList) > 1: static = False
if static:
    print 'All static runs\n'
else:
    for run in runList:
        displayRun(run)      # check if static or not
displayFinal(runList)

fitting = raw_input("Birch fitting? (y/n): ")
# add general minimization option (E vs ayz)
if 'y' in fitting or 'Y' in fitting:
    fins = finalValues(runList)
    energies = fins[1]
    volumes = fins[2]
    fitBirch(energies,volumes)



# also no need to display for unrelaxed runs, check length of energy list??
# print relaxed or unrelaxed
# directory must end in /
# testing
#data = parseResults('/Users/Jonas/Google Drive/Laspa_JLK/computation_files/GSFE/full/')
