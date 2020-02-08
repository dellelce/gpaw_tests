'''
    Multiple executions from modifying json object
'''

import json
from auto import execution, frange
from pprint import pprint as p
import gpaw

inFile = '01.json'

with open(inFile, 'rb') as f:
    j = json.loads(f.read())

p(j)

for _s in frange(2.67, 2.85, 0.01):

    j['xyz'][3] = ['C', _s, 0.0, 0.0]

    x = execution(j)

    p(x.inp)
    p(x.gpts)

    try:
        x.run(str(_s))
        p(x.atoms.get_forces())
    except gpaw.KohnShamConvergenceError:
        pass

## EOF ##
