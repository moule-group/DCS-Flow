import os
from pathlib import Path
from cnss.relax import relax
from cnss.phonons import phonons
from cnss.oclimax import oclimax
from cnss import read_json, write_json, get_default_parameters
from ase.utils.timing import Timer

class CLICommand:
    'Workflow to relax structure, calculate phonons and calculate INS spectrum'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('--params',
            help='JSON file with parameters for workflow',
            default='workflow_params.json')
        add('--get-params',
            action='store_true',
            help='Get JSON file with default parameters for workflow')

    @staticmethod
    def run(args):
        if args.get_params:
            write_params()
            return
        
        dct = {}
        if Path(args.params).is_file():
            dct = read_json(args.params)
        workflow(dct)

def workflow(dct=None):
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
