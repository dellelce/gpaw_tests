
#from ase.lattice.spacegroup import crystal
from ase.spacegroup import crystal
#from parcalc import ClusterVasp, ParCalculate

import ase.units as units
import numpy
#import matplotlib.pyplot as plt

## from auto
from ase.build import bulk
from ase.io import read
from gpaw import GPAW
from gpaw.utilities import h2gpts
from gpaw.mixer import BroydenMixer
from gpaw.mixer import Mixer

from time import time

# for stress calculation - from example:
# https://wiki.fysik.dtu.dk/gpaw/exercises/stress/stress.html

from ase.optimize.bfgs import BFGS
from ase.constraints import UnitCellFilter
from gpaw import PW

a = 4.194
cryst = crystal(['Mg', 'O'],
                [(0, 0, 0), (0.5, 0.5, 0.5)],
                spacegroup=225,
                cellpar=[a, a, a, 90, 90, 90])

print(cryst)

'''
   _calc = GPAW(
                xc=_fnl,
                txt=self.filename(force_id),
                nbands=self.nbands,
                mixer=self.mixer(),
                maxiter=self.max_iterations,
                gpts=self.gpts
               )
'''

calc = GPAW(
            xc='PBE',
            txt='ela_' + str(int(time())) + '.txt'
           )

#

cryst.set_calculator(calc)

# Set the calculation parameters
#calc.set(prec = 'Accurate', xc = 'PBE', lreal = False,
#calc.set(xc = 'PBE', lreal = False,
calc.set(xc = 'PBE', 
         mode=PW(400),
#            nsw=30, ediff=1e-8, ibrion=2, kpts=[3,3,3])
#            ediff=1e-8, ibrion=2, kpts=[3,3,3])
#            ibrion=2, kpts=[3,3,3])
             kpts=[3,3,3])

# Set the calculation mode first.
# Full structure optimization in this case.
# Not all calculators have this type of internal minimizer!
#calc.set(isif=3)

# from ase/GPAW stress example

print ("Potential Energy: " + str(cryst.get_potential_energy()))

uf = UnitCellFilter(cryst)
relax = BFGS(uf)
relax.run(fmax=0.05)

stress = cryst.get_stress()
print(stress)

a=cryst.get_cell()[0,0]

print ("optimized lattice constant: "+a)




### EOF ###
