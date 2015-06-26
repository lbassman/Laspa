# Jonas Kaufman jlkaufman@gmail.com
# June 22 2014

from pylab import *
from Cell import *
import math

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

#===========================================================================
# MAIN PROGRAM
#===========================================================================
"""
element = raw_input('Element: ')
a = float(raw_input('Lattice parameter: '))
faultType = raw_input('Stacking fault type (P/I/E): ')
n = int(raw_input('Number of 3-layer stacks on either side of fault: '))
"""
element = 'Cu'
a = 3.63548
ngap = 3

sequence = 'ABCABCABCABC'
layers = list(sequence)
nlayers = len(layers)

ax = math.sqrt(2)/2
ay = math.sqrt(6)/2
az = math.sqrt(3)/3
bx,by,bz = ax,ay,az*(nlayers+ngap)

cell = Cell() # Cartesian coordinates by default
cell.setHeader(element+' fcc')
cell.setElements([element]) # fix this
cell.setA0(a)
cell.setCoordinateSystem('Direct')
cell.setLatticeVectors([[bx,0,0],[0,by,0],[0,0,bz]])

zpoint = 0
while layers:
    next = layers.pop(0)
    if next == 'A':
        addA(cell,zpoint)
    elif next == 'B':
        addB(cell,zpoint)
    else: # next == 'C'
        addC(cell,zpoint)
    zpoint += (1.0/(nlayers+ngap))

halfz = .5*(nlayers-1)/(nlayers+ngap) # z position right between the surface layers

cell.setSiteMobilities(False,False,True)
cell.sendToPOSCAR()

n_points = 5
displacements = np.linspace(0.0, 1.0, n_points)

# ISF
for disp in displacements:    
    cell = Cell().loadFromPOSCAR()
    for element in cell.sites:
        for i in range(cell.numberOfElements()):
            for j in range(cell.numberOfAtomsOfElement(i)):
                    x,y,z = cell.sites[i][j].position
                    if z > halfz:
                        cell.sites[i][j].move([x, y - disp*(1.0/6.0), z])
    cell.setSiteMobilities(False,False,True) # add conditions for this?
    disp = ('%.3f'%disp).replace('.','')
    cell.sendToPOSCAR('POSCAR_I_%s'%disp)

halfz = .5*(nlayers-3)/(nlayers+ngap)

# ESF
for disp in displacements[1:]:    
    cell = Cell().loadFromPOSCAR('POSCAR_I_1000')
    for element in cell.sites:
        for i in range(cell.numberOfElements()):
            for j in range(cell.numberOfAtomsOfElement(i)):
                    x,y,z = cell.sites[i][j].position
                    if z < halfz:
                        cell.sites[i][j].move([x, y + disp*(1.0/6.0), z])
    cell.setSiteMobilities(False,False,True) # add conditions for this?
    disp = ('%.3f'%disp).replace('.','')
    cell.sendToPOSCAR('POSCAR_E_%s'%disp)
