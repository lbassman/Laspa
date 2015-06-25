# Jonas Kaufman jlkaufman@gmail.com
# June 22 2014
# Script to create POSCARs for several points along the (111)
# generalized stacking fault pathway in an fcc crystal
# and submit all jobs to the queue on Stampede

from Cell import *
import math

#===========================================================================
# MAIN PROGRAM
#===========================================================================
def addA(cell,z):
    cell.addSite(0,Site([0,0,z]))
    cell.addSite(0,Site([0.5,0.5,z]))
def addB(cell,z):
    cell.addSite(0,Site([0,(2.0/6.0),z]))
    cell.addSite(0,Site([0.5,(5.0/6.0),z]))
def addC(cell,z):
    cell.addSite(0,Site([0,(4.0/6.0),z]))
    cell.addSite(0,Site([0.5,(1.0/6.0),z]))

#===========================================================================
# MAIN PROGRAM
#===========================================================================
element = raw_input('Element: ')
a = float(raw_input('Lattice parameter: '))
faultType = raw_input('Stacking fault type (P/I/E): ')
n = int(raw_input('Number of 3-layer stacks on either side of fault: '))

sequence = n*'ABC'
if faultType == 'I':
    sequence += 'BC'
elif faultType == 'E':
    sequence += 'ACBC'
sequence += n*'ABC'
layers = list(sequence)
nlayers = len(layers)
print sequence

ax = math.sqrt(2)/2
ay = math.sqrt(6)/2
az = math.sqrt(3)/3
bx,by,bz = ax,ay,az*nlayers

cell = Cell() # Cartesian coordinates by default
cell.setHeader(element+' fcc')
cell.setElements([element]) # fix this
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
    zpoint += (1.0/nlayers)

cell.setSiteMobilities(False,False,True)
cell.sendToPOSCAR()
