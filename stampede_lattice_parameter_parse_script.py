# Jonas Kaufman jlkaufman@hmc.edu
# March 7, 2015
# Script to parse VASP output files for energy values
# for a set of lattice parameters for cubic structures
# Based on a script by Josh Sanz

import os
import fnmatch

# home and work directories
HOME = '/home1/03022/bassman/Jonas/'
WORK = '/work/03022/bassman/Jonas/'

def getEnergy(file):
    """ parses an OUTCAR file and pulls out the free energy of the system """
    f = open(file,'r')
    while True:
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' in f.readline():
            f.readline()    # line of dashes
            energyLine = f.readline().split()
            energy = float(energyLine[4])
            break
    return energy
    
def getCellSize(file):
    """ returns the volume of a cell parsed from OUTCAR """
    f = open(file,'r')
    while True:
        if 'VOLUME and BASIS-vectors are now :' in f.readline():
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            volumeLine = f.readline().split()
            V = float(volumeLine[4])
            for i in range(6):
                f.readline()    # text
            a = float(f.readline().split()[1]) # perpendicular lattice vectors
            # the lattice vector length in loading direction is [0]
            break
    return (V,a)
    ## check that this is really reading the perpendicular size

def parseResults(directory):
    """ finds the energy of each simulation in a given results directory """
    EList = []
    VList = []
    aList = []
    tList = []
    for dir in os.listdir(directory): # 0.XXXXX_main or 0.XXXXX_c*_0.XXXXX for example
        for file in os.listdir(directory + dir):        
                if fnmatch.fnmatch(file,'OUTCAR'):
                    EList += [getEnergy(directory+dir+'/'+file)]
                    sizeData = getCellSize(directory+dir+'/'+file)
                    VList += [sizeData[0]]
                    aList += [sizeData[1]]
        dirInfo = dir.split('_')
        # eN = dirInfo[0] # main strain amount
        if len(dirInfo) == 2:
            tList += [0.0] # additional strain amount
        else:
            tList += [float(dirInfo[2])]          
    return [EList,VList,aList,tList] 


#===========================================================================
# MAIN PROGRAM
#===========================================================================
RESULTS = HOME # where all the results folders are located

print 'GS results:'
data = parseResults(RESULTS+'GS_results/')
print data
EList = data[0]
VList = data[1]



#===========================================================================
# fitting the Birch-Murnaghan EOS to data
# http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/appendices/appendix-eos.html
#===========================================================================

from pylab import * #this includes numpy as np!
from scipy.optimize import leastsq

# raw data from 2.2.3-al-analyze-eos.py
v = np.array(VList)
e = np.array(EList)

#make a vector to evaluate fits on with a lot of points so it looks smooth
vfit = np.linspace(min(v),max(v),100)

### fit a parabola to the data
# y = ax^2 + bx + c
a,b,c = polyfit(v,e,2) #this is from pylab

'''
the parabola does not fit the data very well, but we can use it to get
some analytical guesses for other parameters.

V0 = minimum energy volume, or where dE/dV=0
E = aV^2 + bV + c
dE/dV = 2aV + b = 0
V0 = -b/2a

E0 is the minimum energy, which is:
E0 = aV0^2 + bV0 + c

B is equal to V0*d^2E/dV^2, which is just 2a*V0

and from experience we know Bprime_0 is usually a small number like 4
'''

#now here are our initial guesses.
v0 = -b/(2*a)
e0 = a*v0**2 + b*v0 + c
b0 = 2*a*v0
bP = 4

#now we have to create the equation of state function
def Murnaghan(parameters,vol):
    '''
    given a vector of parameters and volumes, return a vector of energies.
    equation From PRB 28,5480 (1983)
    '''
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    
    E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)

    return E

# and we define an objective function that will be minimized
def objective(pars,y,x):
    #we will minimize this function
    err =  y - Murnaghan(pars,x)
    return err

x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function

murnpars, ier = leastsq(objective, x0, args=(e,v)) #this is from scipy

#now we make a figure summarizing the results
plot(v,e,'ro')
plot(vfit, a*vfit**2 + b*vfit + c,'--',label='parabolic fit')
plot(vfit, Murnaghan(murnpars,vfit), label='Murnaghan fit')
xlabel('Volume ($\AA^3$)')
ylabel('Energy (eV)')
legend(loc='best')

#add some text to the figure in figure coordinates
ax = gca()
text(0.4,0.5,'Min volume = %1.2f $\AA^3$' % murnpars[3],
     transform = ax.transAxes)
text(0.4,0.4,'Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa' % (murnpars[1],
                                                              murnpars[1]*160.21773)
     , transform = ax.transAxes)
savefig('a-eos.png')
show()

print 'initial guesses  : ',x0
print 'fitted parameters: ', murnpars
   