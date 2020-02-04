#!/usr/bin/env python3

'''
 File:         first.py
 Created:      220517
 Description:  description for first.py
'''

from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton

system = Atoms('H2', positions=[[0.0, 0.0, 0.0],
                                [0.0, 0.0, 1.0]])
calc = EMT()

system.set_calculator(calc)

opt = QuasiNewton(system, trajectory='h2.emt.traj')

opt.run(fmax=0.0001)

## EOF ##
