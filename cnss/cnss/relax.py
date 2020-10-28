import os
from ase.optimize import BFGS
from ase.io import read
from cnss import mkdir, chdir, out, done, isdone

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
        add('--fmax',
            help='Convergence criteria for forces',
            default=0.01,
            type=float)


    @staticmethod
    def run(args):
        relax(args.krelax, args.fmax, args.geo, args.calc)

def relax_structure(krelax, fmax, geo, mode):
    if isdone('relax'):
        return
    else:
        atoms = read(geo)
        formula = atoms.get_chemical_formula()

        if mode == 'dftbp':
            from ase.calculators.dftb import Dftb
            calculator = Dftb(label=formula,
                              atoms=atoms,
                              Driver_='ConjugateGradient',
                              Driver_MovedAtoms='1:-1',
                              Driver_MaxForceComponent=fmax,
                              Driver_MaxSteps=100,
                              # Driver_LatticeOpt='Yes',
                              # Driver_Isotropic='Yes',
                              # Driver_AppendGeometries='Yes',
                              kpts=krelax,
                              Hamiltonian_SCC='Yes',
                              Hamiltonian_MaxAngularMomentum_='',
                              Hamiltonian_MaxAngularMomentum_C='p',
                              Hamiltonian_MaxAngularMomentum_H='s')

        elif mode == 'chimes':
            from ase.calculators.dftb import Dftb
            from cnss.chimes import run_md_input
            folder = os.getcwd()            
            run_md_input(folder + '/..')
            calculator = Dftb(label=formula,
                              atoms=atoms,
                              Driver_='ConjugateGradient',
                              Driver_MovedAtoms='1:-1',
                              Driver_MaxForceComponent=fmax,
                              Driver_MaxSteps=100,
                              kpts=krelax,
                              Hamiltonian_ChIMES='Yes',
                              Hamiltonian_SCC='Yes',
                              Hamiltonian_MaxAngularMomentum_='',
                              Hamiltonian_MaxAngularMomentum_C='p',
                              Hamiltonian_MaxAngularMomentum_H='s')

        elif mode == 'vasp':
            from ase.calculators.vasp import Vasp
            calculator = Vasp(kpts=krelax,
                              prec='Accurate',
                              encut=550,
                              nsw=100,
                              isif=2,
                              ismear=0,
                              ibrion=1,
                              ediff=1e-8,
                              ediffg=-fmax,
                              sigma=0.1,
                              nwrite=1,
                              ncore=16,
                              lreal=False,
                              lcharg=False,
                              lwave=False,
                              xc='pbe',
                              gamma=True)
        else:
            raise NotImplementedError('{} calculator not implemented' .format(mode))
            
        atoms.set_calculator(calculator)
        atoms.get_potential_energy()
        done('relax')

def find_geo(folder):
    import glob
    
    geo = glob.glob(folder + '/*.cif') + \
          glob.glob(folder + '/*.gen') + \
          glob.glob(folder + '/*.sdf') + \
          glob.glob(folder + '/*.xyz')
    geo = geo[0]

    return geo
        
def relax(krelax=[6, 6, 6], fmax=0.01, geo=None, calc='dftbp'):
    folder = os.getcwd()

    if geo:
        try:
            geo = folder + '/../' + geo
        except:
            geo = folder + '/' + geo
    else:
        geo = find_geo(folder)
        
    mkdir(folder + '/1-optimization')
    with chdir(folder + '/1-optimization'):
        with out('relax'):
            relax_structure(krelax=krelax, fmax=fmax, geo=geo, mode=calc)

if __name__ == '__main__':
    relax()
