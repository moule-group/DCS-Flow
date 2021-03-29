import os
from pathlib import Path
from dcs.relax import relax
from dcs.phonons import phonons
from dcs.oclimax import oclimax
from dcs import read_json, write_json, get_default_parameters
from ase.utils.timing import Timer

class CLICommand:
    'Workflow to relax structure, calculate phonons and calculate INS spectrum'

    @staticmethod
    def add_arguments(parser):
        """Sets up command line to run workflow script i.e. recognize arguments and commands.

        Args:
            parser (argparse): Arguments to be added. 
        """
        add = parser.add_argument
        add('--params',
            help='JSON file with parameters for workflow',
            default='workflow_params.json')
        add('--get-params',
            action='store_true',
            help='Get JSON file with default parameters for workflow')

    @staticmethod
    def run(args):
        """Runs workflow function using command line arguments. 

        Args:
            args (argparse): Command line arguments added to parser using the function add_arguments.
        """
        if args.get_params:
            write_params()
            return
        
        dct = {}
        if Path(args.params).is_file():
            dct = read_json(args.params)
        workflow(dct)

def workflow(dct=None):
    """Calls all workflow functions (relax, phonon, oclimax) with a timer using specified parameters in workflow_params.json, else with default parameters.

    Args:
        dct (dict, optional): Specified parameters for relax, phonons, and oclimax functions. Defaults to None.
    """
    timer = Timer()
    
    if dct:
        with timer('relaxation'):
            relax(dct['krelax'], dct['fmax'], dct['geo'], dct['calc'])
        timer.write()
        with timer('phonon calculation'):
            phonons(dct['dim'], dct['kforce'], dct['mesh'], dct['calc'])
        timer.write()
        with timer('oclimax calculation'):
            oclimax(dct['params'], dct['task'], dct['e_unit'])
        timer.write()

    else:
        with timer('relaxation'):
            relax()
        timer.write()
        with timer('phonon calculation'):
            phonons()
        timer.write()
        with timer('oclimax calculation'):
            oclimax()
        timer.write()

def write_params():
    """Writes json file with default arguments from relax, phonons, and oclimax scripts. 
    """
    relax_params = get_default_parameters(relax)
    phonons_params = get_default_parameters(phonons)
    oclimax_params = get_default_parameters(oclimax)

    params = {**relax_params, **phonons_params, **oclimax_params}
    write_json('workflow_params.json', params)

if __name__ == '__main__':
    try:
        dct = read_json('workflow_params.json')
        workflow(dct)
    except:
        workflow()
