import os
import glob
from ase.optimize import BFGS
from ase.io import read
from cnss import mkdir, chdir, out


class CLICommand:
    'Optimize structure'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
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

def relax_done(fmax):
    import pandas as pd
    from pandas.io import common

    try:
        df = pd.read_csv('relax.out', sep='\s+')
        if (df['fmax'] < fmax).any():
            return True
        else:
            return False
    except:
        return False

    
def relax_structure(krelax, fmax, geo, calc):
    if relax_done(fmax = fmax):
        return
    else:
        atoms = read(geo)
        formula = atoms.get_chemical_formula()

        if calc == 'dftbp':
            from ase.calculators.dftb import Dftb
            calc = Dftb(label=formula,
                        atoms=atoms,
                        kpts=krelax,
                        Hamiltonian_SCC='Yes',
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_H='s'))
        elif calc == 'vasp':
            from ase.calculators.dftb import Vasp
            calc = Vasp(kpts=krelax,
                        encut=520,
                        prec='Accurate',
                        nwrite=1,
                        lcharg=False,
                        lwave=False,
                        xc='optPBE-vdw',
                        gamma=True)
        else:
            raise NotImplementedError('{} calculator not implemented' .format(calc))
            
        atoms.set_calculator(calc)
        opt = BFGS(atoms, trajectory= formula + '.traj')
        opt.run(fmax=fmax)

        
def relax(krelax=[6, 6, 6], fmax=0.01, geo=None, calc='dftbp'):
    folder = os.getcwd()
    if not geo:
        geo = glob.glob(folder + '/*.cif') + \
            glob.glob(folder + '/*.gen') + \
            glob.glob(folder + '/*.sdf') + \
            glob.glob(folder + '/*.xyz')

    mkdir(folder + '/1-optimization')
    with chdir(folder + '/1-optimization'):
        with out('relax'):
            relax_structure(krelax=krelax, fmax=fmax, geo=geo[0], calc=calc)

if __name__ == '__main__':
    relax()
