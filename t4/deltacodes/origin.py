'''
from pdf
GPAW09(g)/GPAW
GPAW 0.9 PAW dataset / GPAW 0.10.0 grid-based
script for the Atomic Simulation Environment (ASE)
S. R. Bahn, K. W. Jacobsen, Comput. Sci. Eng. 4, 56{66 (2002)
https://wiki.fysik.dtu.dk/ase/
'''

import os
import sys
import time
import numpy as np
import ase.db
from ase.units import Rydberg
from ase.utils import opencew
from ase.calculators.calculator import kpts2mp
from ase.io.trajectory import PickleTrajectory
from ase.test.tasks.dcdft import DeltaCodesDFTCollection as Collection
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.utilities import h2gpts

collection = Collection()

if len(sys.argv) == 1:
 names = collection.names
else:
 names = [sys.argv[1]]

c = ase.db.connect('dcdft_gpaw_fd.db')

mode = 'fd'
e = 0.08 # h -> gpts
kptdensity = 16.0 # this is converged
width = 0.01
relativistic = True
constant_basis = True

if relativistic:
linspace = (0.98, 1.02, 5) # eos numpy's linspace
else:
linspace = (0.92, 1.08, 7) # eos numpy's linspace
linspacestr = ''.join([str(t) + 'x' for t in linspace])[:-1]
code = 'gpaw' + '-' + mode + str(e) + '_c' + str(constant_basis) + '_e' + linspacestr
code = code + '_k' + str(kptdensity) + '_w' + str(width)
code = code + '_r' + str(relativistic)
for name in names:
# save all steps in one traj file in addition to the database
# we should only used the database c.reserve, but here
# traj file is used as another lock ...



