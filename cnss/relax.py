import os
from ase.optimize import BFGS
from ase.io import read
from cnss import mkdir, chdir, out, done, isdone

class CLICommand:
    """Sets up command line to run relax. 
    """

    @staticmethod
    def add_arguments(parser):
        """Sets up command line to run relax script i.e. recognize arguments and commands.

        Args:
            parser (argparse): Arguments to be added. 
        """
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
        """Runs relax.py functions using command line arguments. 

        Args:
            args (argparse): Command line arguments added to parser using the function add_arguments.
        """
        relax(args.krelax, args.fmax, args.geo, args.calc)

def relax_structure(krelax, fmax, geo, mode):
    """Imports specified calculator and calculates optimized structure. [q] or just relaxes structure

    Args:
        krelax (list): Number of k points for relaxation. 
        fmax (float): Maximum allowed force for convergence between atoms.
        geo (str): Geometry file or structure. Allowed file types are cif, gen, sdf, or xyz. 
        mode (str): Calculator used. Options are 'dftbp', 'chimes', or 'vasp'. 

    Raises:
        NotImplementedError: If specified calculator/mode is not one of the above three, error raised. 
    """
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
                              Driver_empty='MaxForceComponent[eV/AA] = {}' .format(fmax),
                              Driver_MaxSteps=1000,
                              Driver_LatticeOpt='Yes',
                              Driver_Isotropic='Yes',
                              Driver_AppendGeometries='Yes',
                              kpts=krelax,
                              Hamiltonian_SCC='Yes',
                              Hamiltonian_MaxAngularMomentum_='',
                              Hamiltonian_MaxAngularMomentum_C='p',
                              Hamiltonian_MaxAngularMomentum_O='p',
                              Hamiltonian_MaxAngularMomentum_H='s',
                              Hamiltonian_MaxAngularMomentum_N='p',
                              Hamiltonian_MaxAngularMomentum_S='d')

        elif mode == 'chimes':
            from ase.calculators.dftb import Dftb
            from cnss.chimes import run_md_input
            folder = os.getcwd()            
            run_md_input(folder + '/..')
            calculator = Dftb(label=formula,
                              atoms=atoms,
                              Driver_='ConjugateGradient',
                              Driver_MovedAtoms='1:-1',
                              Driver_empty='MaxForceComponent[eV/AA] = {}' .format(fmax),
                              Driver_MaxSteps=1000,
                              Driver_LatticeOpt='Yes',
                              Driver_Isotropic='Yes',
                              Driver_AppendGeometries='Yes',
                              kpts=krelax,
                              Hamiltonian_ChIMES='Yes',
                              Hamiltonian_SCC='Yes',
                              Hamiltonian_MaxAngularMomentum_='',
                              Hamiltonian_MaxAngularMomentum_C='p',
                              Hamiltonian_MaxAngularMomentum_O='p',
                              Hamiltonian_MaxAngularMomentum_H='s',
                              Hamiltonian_MaxAngularMomentum_N='p',
                              Hamiltonian_MaxAngularMomentum_S='d')

        elif mode == 'vasp':
            from ase.calculators.vasp import Vasp
            calculator = Vasp(kpts=krelax,
                              prec='Accurate',
                              encut=520,
                              nsw=100,
                              isif=3,
                              ismear=0,
                              ibrion=1,
                              ediff=1e-8,
                              ediffg=-fmax,
                              sigma=0.1,
                              nwrite=1,
                              ncore=16,
                              lreal='Auto',
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
    """Searches the current working directory for the molecular structure file. 
        Allowed file types are cif, gen, sdf, or xyz. 

    Args:
        folder (str): The current working directory.

    Returns:
        str: Molecular structure file. 
    """
    import glob
    
    geo = glob.glob(folder + '/*.cif') + \
          glob.glob(folder + '/*.gen') + \
          glob.glob(folder + '/*.sdf') + \
          glob.glob(folder + '/*.xyz')
    geo = geo[0]

    return geo
        
def relax(krelax=[6, 6, 6], fmax=0.01, geo=None, calc='dftbp'):
    """Finds geo file and runs relax_structure funtion using specified calculator. 

    Args:
        krelax (list, optional): Number of k points for relaxation. Defaults to [6, 6, 6].
        fmax (float, optional): Maximum allowed force for convergence between atoms. Defaults to 0.01.
        geo (str, optional): Geometry file or structure. 
            Allowed file types are cif, gen, sdf, or xyz. Defaults to None.
        calc (str, optional): Calculator used. Options are 'dftbp', 'chimes', or 'vasp'. Defaults to 'dftbp'.
    """
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
