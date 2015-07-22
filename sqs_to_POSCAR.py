#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  July 20, 2015
#  Script to convert a bestsqs.out file from ATAT to a POSCAR file for use  
#  with VASP
#  Uses direct coordinates
#==============================================================================
"""
Make script executable using 'chmod +x _____.py' to call as bash script
Requires Cell.py
"""
from Cell import *

cell = Cell().loadFromSQS('bestsqs.out')
cell.sendToPOSCAR('POSCAR')