
'''
  from:
  https://wiki.fysik.dtu.dk/gpaw/tutorials/gw_tutorial/gw_tutorial.html
'''

from ase.build import bulk
from ase.io import read
from gpaw import GPAW
from gpaw.utilities import h2gpts

import json

class execution(object):
 
 def __init__(self, name):
  pass



def main():
 atoms = read('g1_001.xyz')
 a = 8.467
 atoms.cell=(a-0.5, a, 3)
 atoms.pbc=False
 atoms.center()
 
 print(atoms.get_cell())
 
 calc = GPAW(
            xc='M06-L',
            txt='C_groundstate_M06L_9.txt',
            nbands=60,
            gpts=h2gpts(0.20,atoms.get_cell(),idiv=8)
           )

 atoms.set_calculator(calc)
 atoms.get_potential_energy()

 #calc.diagonalize_full_hamiltonian()       # determine all bands
 calc.write('C_groundstate.gpw','all')    # write out wavefunctions


## call main

main()

## EOF ##
