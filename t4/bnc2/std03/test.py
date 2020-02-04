'''
    Multiple executions from modifying json object
'''

import json
from auto import execution, frange
from pprint import pprint as p
import gpaw
import os


'''
     [ [ "N",  0.0000, 0.0000, 0.0000 ],
                 [ "C",  1.2500, 0.0000, 0.0000 ],
                 [ "C",  2.6500, 0.0000, 0.0000 ]
'''

inFile = '01.noboron.fire.pbe.json'
os.environ['OMP_NUM_THREADS'] = '1'

with open(inFile,'rb') as f:
 j = json.loads(f.read())

p(j)
step=0.2

for _x0 in frange(-0.14,0.15,step):
 j['xyz'][0] = ['N',_x0,0.0,0.0]

 for _x1 in frange(1.10, 1.40, step):
  j['xyz'][1] = ['C',_x1,0.0,0.0]

  for _x2 in frange(2.40, 2.80, step):
   j['xyz'][2] = ['C',_x2,0.0,0.0]

   print("====")
   p(j['xyz'])
   x = execution(j)
   p(x.inp)
   p(x.gpts)

   label=str(_x0)+'_'+str(_x1)+'_'+str(_x2)
   label.replace('.','')

   try:
    x.run(label)
    p(x.atoms.get_forces())
   except gpaw.KohnShamConvergenceError:
    print('*NOT CONVERGED*')
   finally:
    print("   Elapsed => " + str(x.lastElapsed))

## EOF ##
