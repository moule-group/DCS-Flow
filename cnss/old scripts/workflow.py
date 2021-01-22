import os
from pathlib import Path
from cnss.relax import relax
from cnss.phonons import phonons
from cnss.oclimax import oclimax
from cnss import read_json, write_json, get_default_parameters
from ase.utils.timing import Timer

class CLICommand:
    # class to define help function
    'Workflow to relax structure, calculate phonons and calculate INS spectrum'

    @staticmethod # decorator
    def add_arguments(parser):
        # set up arguments and parse together
        add = parser.add_argument
        add('--params',
            help='JSON file with parameters for workflow',
            default='workflow_params.json')
        add('--get-params',
            action='store_true',
            help='Get JSON file with default parameters for workflow')

    @staticmethod
    def run(args):
        # executes the command, checks for specified params and writes them to the json file
        if args.get_params:
            write_params()
            return
        
        dct = {}
        if Path(args.params).is_file():
        # if any json file exists in the folder with the below specified name, will set dct equal to it
            dct = read_json(args.params)
        workflow(dct)

def workflow(dct=None): # default value of dct=none (evaluated as false), not completely clear on value of dct
    # 
    timer = Timer()
    
    if dct: # if workflow(dct) run with dct argument (what should dct be? params file)
        with timer('relaxation'):
            relax(dct['krelax'], dct['fmax'], dct['geo'], dct['calc'])
        with timer('phonon calculation'):
            phonons(dct['dim'], dct['kforce'], dct['mesh'], dct['calc'])
        with timer('oclimax calculation'):
            oclimax(dct['params'], dct['task'], dct['e_unit'])
    else: # if there's no json file
        with timer('relaxation'):
            relax()
        with timer('phonon calculation'):
            phonons()
        with timer('oclimax calculation'):
            oclimax()

    timer.write()

def write_params():
    # gets defaults defined in other scripts
    relax_params = get_default_parameters(relax)
    phonons_params = get_default_parameters(phonons)
    oclimax_params = get_default_parameters(oclimax)

    # sets params as defaults defined above
    params = {**relax_params, **phonons_params, **oclimax_params}
    # writes params to json file as dict
    write_json('workflow_params.json', params)

if __name__ == '__main__':
    try:
        # tries to reads params if the user provided them, if not then it'll go to exceot (won't return error though)
        dct = read_json('workflow_params.json')
        # runs workflow using json file (can edit and will use those arguments)
        workflow(dct)
    except: # what will run if no json file exists (default arguments) provided by functions above
        workflow()
