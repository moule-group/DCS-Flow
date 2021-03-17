import os
from ase.optimize import BFGS
from ase.io import read
from cnss import mkdir, chdir, out, done, isdone

class CLICommand:
    'Optimize structure'

    @staticmethod
    def add_arguments(parser):
        """Sets up command line to run relax script i.e. recognize arguments and commands.

        Args:
            parser (argparse): Arguments to be added. 
        """
        add = parser.add_argument
        add('--calc', help='Calculator used. Options are dftbp, vasp or castep', default='dftbp')
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
    """Defines arguments for specified calculator and optimizes the structure.
    
    Args:
        krelax (list): Number of k points for relaxation. 
        fmax (float): Maximum allowed force for convergence between atoms.
        geo (str): Geometry file or structure. Allowed file types are .cif, .gen, .sdf, or .xyz. 
        mode (str): Calculator used. Options are 'dftbp', 'chimes', 'vasp', or 'castep'. 

    Raises:
        NotImplementedError: If specified calculator is not an option for mode, error raised. 
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
                              Driver_='LBFGS',
                              Driver_MovedAtoms='1:-1',
                              Driver_empty='MaxForceComponent[eV/AA] = {}' .format(fmax),
                              Driver_MaxSteps=1000,
                              Driver_LatticeOpt='Yes',
                              Driver_Isotropic='Yes',
                              kpts=krelax,
                              Hamiltonian_SCC='Yes',
                              Hamiltonian_SCCTolerance=1e-7,
                              Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}' .format(T=5),
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
                              Driver_='LBFGS',
                              Driver_MovedAtoms='1:-1',
                              Driver_empty='MaxForceComponent[eV/AA] = {}' .format(fmax),
                              Driver_MaxSteps=1000,
                              Driver_LatticeOpt='Yes',
                              Driver_Isotropic='Yes',
                              kpts=krelax,
                              Hamiltonian_ChIMES='Yes',
                              Hamiltonian_SCC='Yes',
                              Hamiltonian_SCCTolerance=1e-7,
                              Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}' .format(T=5),
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
                              nsw=1000,
                              isif=3,
                              ismear=0,
                              ibrion=1,
                              ediff=1e-8,
                              ediffg=-fmax,
                              sigma=0.05,
                              nwrite=1,
                              ncore=16,
                              lreal=False,
                              lcharg=False,
                              lwave=False,
                              xc='pbe',
                              gamma=True)

        elif mode == 'castep':
            import ase.calculators.castep

            calculator = ase.calculators.castep.Castep()
            directory = '../1-optimization'
            calculator._export_settings = True
            calculator._directory = directory
            calculator._rename_existing_dir = False
            calculator._export_settings = True
            calculator._label = 'relax'
            calculator._set_atoms = True
            calculator.param.task = 'GeometryOptimization'
            calculator.param.xc_functional = 'PBE'
            calculator.param.basis_precision = 'MEDIUM'
            calculator.param.geom_method = 'BFGS'
            calculator.param.cut_off_energy = 520
            calculator.param.num_dump_cycles = 0
            calculator.param.geom_force_tol = fmax
            calculator.param.geom_energy_tol = 1e-5
            calculator.param.geom_disp_tol = 1e-3
            calculator.param.elec_energy_tol = 1e-8
            calculator.param.geom_max_iter = 1000
            calculator.cell.kpoint_mp_grid = krelax            
            
        else:
            raise NotImplementedError('{} calculator not implemented' .format(mode))
            
        atoms.set_calculator(calculator)
        atoms.get_potential_energy()
        done('relax')

def find_geo(folder):
    """Searches the current working directory for the molecular structure file. 
        Allowed file types are .cif, .gen, .sdf, or .xyz. 

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
        
def relax(krelax=[6, 6, 6], fmax=0.05, geo=None, calc='dftbp'):
    """Finds the geometry file and optimizes the structure using the specified calculator 
    (Populates 1-optimization folder with results). 

    Args:
        krelax (list, optional): Number of k points for relaxation. Defaults to [6, 6, 6].
        fmax (float, optional): Maximum allowed force for convergence between atoms. Defaults to 0.01.
        geo (str, optional): Geometry file or structure. 
            Allowed file types are .cif, .gen, .sdf, or .xyz. Defaults to None.
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
