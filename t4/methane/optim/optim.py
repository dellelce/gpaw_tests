from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.io import read

system = read('octo_0.xyz')
calc = EMT()

system.set_calculator(calc)

opt = QuasiNewton(system, trajectory='h2.emt.traj')

opt.run(fmax=0.05)
