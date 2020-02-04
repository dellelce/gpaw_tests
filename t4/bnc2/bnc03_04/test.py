'''
    Multiple executions from modifying json object
'''

import json
from auto import execution, frange
from pprint import pprint as p
import gpaw
import os

inFile = '01.json'
os.environ['OMP_NUM_THREADS'] = '1'

with open(inFile,'rb') as f:
 j = json.loads(f.read())

p(j)

for _s in frange(2.86,3,0.01) :

 j['xyz'][3] = ['C',_s,0.0,0.0]

 x = execution(j)

 p(x.inp)
 p(x.gpts)

 try:
  x.run(str(_s))
  p(x.atoms.get_forces())
 except gpaw.KohnShamConvergenceError:
  print('*NOT CONVERGED*')
 finally:
  print("   Elapsed => " + str(x.lastElapsed))

## EOF ##
