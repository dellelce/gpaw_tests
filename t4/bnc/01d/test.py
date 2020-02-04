'''
    Multiple execution from modifying json object
'''

import json
from auto import execution
from pprint import pprint as p

inFile = '00.json'

with open(inFile,'rb') as f:
 j = json.loads(f.read())

for _s in [ 6.6, 6.7, 6.8, 6.9, 7, 7.1 ]: 

 j['cell']['x'] = _s 

 x = execution(j)
 p(x.inp)
 p(x.gpts)

 x.run(str(_s))

 p(x.atoms.get_forces())

## EOF ##
