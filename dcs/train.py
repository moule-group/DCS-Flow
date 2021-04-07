import os
from pathlib import Path
from dcs.relax import relax, find_geo
from dcs.md import md
from dcs.chimes import chimes
from dcs import (read_json, write_json, get_default_parameters,
                  mkdir, chdir)
from ase.utils.timing import Timer
from shutil import copy

class CLICommand:
    'Workflow to run DFT-MD and train ChIMES model'

    @staticmethod
    def add_arguments(parser):
        """Sets up command line to run train script i.e. recognize arguments and commands.

        Args:
            parser (argparse): Arguments to be added. 
        """
        add = parser.add_argument
        add('--params',
            help='JSON file with parameters for training workflow',
            default='train_params.json')
        add('--get-params',
            action='store_true',
            help='Get JSON file with default parameters for training workflow')

    @staticmethod
    def run(args):
        """Runs train function using command line arguments. 

        Args:
            args (argparse): Command line arguments added to parser using the function add_arguments.
        """
        if args.get_params:
            write_params()
            return
        
        dct = {}
        if Path(args.params).is_file():
            dct = read_json(args.params)
        train(dct)

def train(dct=None):
    """Calls functions related to the training workflow (relax, md, chimes) with a timer using specified parameters in train_params.json, else with default parameters. Creates 0-train directory.

    Args:
        dct (dict, optional): JSON file with specified parameters for relax, md, and chimes functions. Defaults to 'train_params.json'.
    """
    timer = Timer()
    folder = os.getcwd()
        
    mkdir(folder + '/0-train')
    geo = find_geo(folder)
    copy(geo, folder + '/0-train/.')
    with chdir(folder + '/0-train'):
        if dct:
            with timer('relaxation'):
                relax(dct['krelax'], dct['fmax'], dct['geo'], dct['calc'])
            timer.write()
            with timer('molecular dynamics'):
                md(dct['optgeo'], dct['calc'], dct['T'], dct['md_size'],
                   dct['steps'], dct['time_step'], dct['dump_interval'])
            timer.write()
            with timer('force matching'):
                chimes(dct['trajfile'], dct['b2'], dct['b3'], dct['T'])
            timer.write()

        else:
            with timer('relaxation'):
                relax(calc='vasp')
            timer.write()
            with timer('molecular dynamics'):
                md()
            timer.write()
            with timer('force matching'):
                chimes()
            timer.write()

    copy(folder + '/0-train/3-chimes/params.txt', folder + '/params.txt')

def write_params():
    """Writes a json file with the arguments and default values for realx, md, and chimes functions.
    """
    relax_params = get_default_parameters(relax)
    md_params = get_default_parameters(md)
    chimes_params = get_default_parameters(chimes)

    params = {**relax_params, **md_params, **chimes_params}
    write_json('train_params.json', params)

if __name__ == '__main__':
    try:
        dct = read_json('train_params.json')
        train(dct)
    except:
        train()