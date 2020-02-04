
'''
  from:
  https://wiki.fysik.dtu.dk/gpaw/tutorials/gw_tutorial/gw_tutorial.html
'''

from ase.build import bulk
from ase.io import read
from gpaw import GPAW
from gpaw.utilities import h2gpts


def main():
 a = 8.267
 atoms = read('g1_001.xyz')
 atoms.cell=(a, a, a)
 atoms.pbc=False
 atoms.center()
 
 print(atoms.get_cell())
 
 calc = GPAW(
            xc='HYB_MGGA_XC_M06_2X',
            txt='C_groundstate_M062X_HYB_0.txt',
            gpts=h2gpts(0.20,atoms.get_cell(),idiv=8)
           )

 atoms.set_calculator(calc)
 atoms.get_potential_energy()

 #calc.diagonalize_full_hamiltonian()       # determine all bands
 calc.write('C_groundstate.gpw','all')    # write out wavefunctions


## call main

main()

## EOF ##
