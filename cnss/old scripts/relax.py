import os
import glob
from ase.optimize import BFGS
from ase.io import read
from cnss import mkdir, chdir, out
import pandas as pd


class CLICommand:
    # adding the argument/subcommand and help explanation
    'Optimize structure'

    @staticmethod
    def add_arguments(parser):
        # defining function add_arguments with argument=parser, not sure where parser comes from
        add = parser.add_argument # parser=object with add_argument attribute
        add('--calc', help='Calculator used. Options are dftbp or vasp', default='dftbp')
        add('--geo', help='Name of geometry file for structure (cif or gen extensions)')
        add('--krelax',
            help='Number of k points for relaxation, e.g., 6 6 6',
            default=[6, 6, 6],
            nargs=3,
            type=int)
        add('--fmax', help='Convergence criteria for forces', default=0.01)


    @staticmethod
    def run(args):
        relax(args.krelax, args.fmax, args.geo, args.calc)

def relax_done(fmax:int): -> bool
    """checks for relax.out file, if relax.out files exists = True """
    try:
        df = pd.read_csv('relax.out', sep='\s+')
        if (df['fmax'] < fmax).any():
            return True
        else:
            return False
    except:
        return False

    
def relax_structure(krelax, fmax, geo, mode): # -> geom optimized struct
    """Uses specified method to optimize molecule geometry """
    if relax_done(fmax=fmax):
        return
    else:
        atoms = read(geo)
        formula = atoms.get_chemical_formula()

        if mode == 'dftbp':
            from ase.calculators.dftb import Dftb
            calculator = Dftb(label=formula,
                              atoms=atoms,
                              kpts=krelax,
                              Hamiltonian_SCC='Yes',
                              Hamiltonian_MaxAngularMomentum_='',
                              Hamiltonian_MaxAngularMomentum_C='p',
                              Hamiltonian_MaxAngularMomentum_H='s')
        elif mode == 'vasp':
            from ase.calculators.vasp import Vasp
            calculator = Vasp(kpts=krelax,
                              encut=520,
                              prec='Accurate',
                              nwrite=1,
                              ncore=16,
                              lreal='Auto',
                              lcharg=False,
                              lwave=False,
                              xc='optpbe-vdw',
                              gamma=True)
        else:
            raise NotImplementedError('{} calculator not implemented' .format(mode))
            
        atoms.set_calculator(calculator)
        opt = BFGS(atoms, trajectory= formula + '.traj')
        opt.run(fmax=fmax)

        
def relax(krelax=[6, 6, 6], fmax=0.01:int, geo=None, calc='dftbp':str): #sets default values, i added the :int and :str but dunno if necessary for this sort of thing
    """ find molecule file, checks in geo opt done, """
    folder = os.getcwd() #sets folder = current wd
    if not geo: # if geo file doesn't exist?, finds molecule file (cif, gen, etc)
        geo = glob.glob(folder + '/*.cif') + \
            glob.glob(folder + '/*.gen') + \
            glob.glob(folder + '/*.sdf') + \
            glob.glob(folder + '/*.xyz')

    mkdir(folder + '/1-optimization') # creates 1-optimization folder
    with chdir(folder + '/1-optimization'):
        with out('relax'): # runs relax_structure
            relax_structure(krelax=krelax, fmax=fmax, geo=geo[0], mode=calc)

if __name__ == '__main__':
    relax()
