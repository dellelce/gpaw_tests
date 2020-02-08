'''
  from:
  https://wiki.fysik.dtu.dk/gpaw/tutorials/gw_tutorial/gw_tutorial.html
'''

from ase.build import bulk
from ase.io import read
from gpaw import GPAW
from gpaw.utilities import h2gpts
from gpaw.mixer import BroydenMixer

#Spacegroup
import spglib

import json
from pprint import pprint

#id generation
import time

### CLASSES ###


class execution(object):
    def __init__(self, executionInput, xyz=None, stateFile='state.json'):
        '''
  Setup a new execution

  Name:  input file
  xyz:   xyz coordinates (for future automated execution)
  stateFile: state file to use for tracking of executions

  TODO: 

  - only ground state in non-periodic cell supported NOW, more to add!
  - inputs from a dict to better handle automated executions
  '''
        if isinstance(executionInput, str):
            self.inp = self.json(executionInput)
            self.inputFileName = executionInput

        if isinstance(executionInput, dict):
            self.inp = executionInput
            self.inputFileName = 'dict'  # TODO: generate meaningful name here (how?)

        self.functional = self.inp.get('functional')

        if self.functional is None:
            self.functional = 'LDA'

        self.atoms = read(self.inp['xyz'])

        self.atoms.cell = (self.inp['cell']['x'], self.inp['cell']['y'],
                           self.inp['cell']['z'])

        # Custom ID: allow to trace parameter set using a custom id
        self.custom_id = self.inp.get('custom_id')

        # Periodic Boundary Conditions
        self.atoms.pbc = False
        self.atoms.center()

        self.nbands = self.inp.get('nbands')

        self.gpts = None

        # Execution state - keep track of execution in a json file
        self.stateFile = stateFile
        self.stateFileLoaded = None  # Content of state file loaded as a dict
        self.state = None

        # Name: to be set by __str__ first time it is run
        self.name = None

        #
        self.logExtension = 'txt'
        self.dataExtension = 'gpw'

    def setup_state(self):
        '''setup state file

     Each Execution State entry had:

     - Generated name
     - Last execution Parameters
     - Original Position
     - Full log name
     - Execution sequence
     - Last Start Time
     - Last End Time
     - Last Elapsed
     - Execution state: Prepared, Started, Completed
     - Convergence: TBD
     - Results: Energy: TBD
  '''

        self.state = {}
        name = self.__str__()
        # incomplete

    #
    def setup_gpts(self):
        self.gpts = h2gpts(self.inp.get('spacing') or 0.20,
                           self.atoms.get_cell(),
                           idiv=self.inp.get('idiv') or 8)

    def filename(self, force_id=None, extension='txt'):
        '''to correctly support this stateFile support must be completed
     as it isn't just add "_0" to __str__() for now
  '''
        base = self.__str__()

        if self.custom_id is not None:
            base = base + '_' + self.custom_id

        if force_id is None:
            _fn = base + '_0.' + extension
        else:
            _fn = '{}_{}.{}'.format(base, force_id, extension)

        return _fn

    #
    @property
    def xyz(self):
        return self.inp['xyz']

    @property
    def spacegroup(self):
        '''international number of determined spacegroup'''
        return spglib.get_symmetry_dataset(spglib.refine_cell(
            self.atoms))['international']

    #
    def __str__(self):
        '''get execution name - to be used for log filename

  Components:

  - Number atoms
  - Volume
  - Spacegroup number
  - Functional name

  '''

        if self.name is None:
            self.name = '{}_{}_{}_{}'.format(
                len(self.atoms), int(self.atoms.get_volume() * 100),
                self.spacegroup, self.functional)

        return self.name

    #
    def pretty(self):
        '''pretty print original parameters'''
        pprint(self.inp)

    #
    def json(self, filename):
        '''load parameters from json'''

        with open(filename, 'rb') as f:
            s = f.read()
            j = json.loads(s)

        return j

    #
    def calc(self, force_id=None):
        '''prepare and launch execution
     force_id: use this id for filename generation
  '''

        self.setup_gpts()

        # generate filenames
        _fn = self.filename(force_id)
        _dn = self.filename(force_id, extension=self.dataExtension)

        print('Calculation with {} log and Input: {}.'.format(
            _fn, self.inputFileName))

        if self.nbands is None:
            _calc = GPAW(xc=self.functional,
                         txt=self.filename(force_id),
                         mixer=BroydenMixer(),
                         gpts=self.gpts)
        else:
            _calc = GPAW(xc=self.functional,
                         txt=self.filename(force_id),
                         nbands=self.nbands,
                         mixer=BroydenMixer(),
                         gpts=self.gpts)

        self.atoms.set_calculator(_calc)
        self.atoms.get_potential_energy()

        _calc.diagonalize_full_hamiltonian()  # determine all bands

        # write out wavefunctions
        _calc.write(_dn, 'all')


###
def main(argv):

    try:
        inputFile = argv[1]
    except:
        inputFile = 'input.json'

    inp = execution(inputFile)

    print(inp)
    inp.calc(str(int(time.time())))


### ORIGINAL ###
def OLDmain():
    atoms = read('g1_001.xyz')
    a = 8.467
    atoms.cell = (a - 0.5, a, 3)
    atoms.pbc = False
    atoms.center()

    print(atoms.get_cell())

    calc = GPAW(xc='M06-L',
                txt='C_groundstate_M06L_9.txt',
                nbands=60,
                gpts=h2gpts(0.20, atoms.get_cell(), idiv=8))

    atoms.set_calculator(calc)
    atoms.get_potential_energy()

    #calc.diagonalize_full_hamiltonian()       # determine all bands
    calc.write('C_groundstate.gpw', 'all')  # write out wavefunctions


## call main

if __name__ == '__main__':
    import sys
    main(sys.argv)

## EOF ##
