import os
import numpy as np
from shutil import copyfile, move
from ase.io import read, write
from dcs import mkdir, chdir, out, done, isdone
from phonopy import Phonopy

class CLICommand:
    'Calculate phonons'

    @staticmethod
    def add_arguments(parser):
        """Sets up command line to run phonopy i.e. recognize arguments and commands.

        Args:
            parser (argparse): Arguments to be added. 
        """
        add = parser.add_argument
        add('--calc',
            help='Calculator used. Options are dftbp, vasp or castep',
            default='dftbp')
        add('--dim',
            help='Dimension of supercell, e. g., 4 4 4',
            default=[4, 4, 4],
            nargs=3,
            type=int)
        add('--kforce',
            help='Number of k points for force calculation, e. g., 1 1 1',
            default=[1, 1, 1],
            nargs=3,
            type=int)
        add('--mesh',
            help='Size of mesh along each dimension, e. g., 8 8 8',
            default=[8, 8, 8],
            nargs=3,
            type=int)


    @staticmethod
    def run(args):
        """Runs phonons function using command line arguments. 

        Args:
            args (argparse): Command line arguments added to parser using the function add_arguments.
        """
        phonons(args.dim, args.kforce, args.mesh, args.calc)


def generate_supercell(dim, mode):
    """Reads in the molecule structure file and creates displacements.

    Args:
        dim (list): Dimensions of the supercell. 
        mode (str): Calculator used for the task. Options are 'dftbp' or 'chimes'. 
    """
    from phonopy.cui.collect_cell_info import collect_cell_info
    from phonopy.interface.calculator import (write_supercells_with_displacements,
                                              get_default_physical_units,
                                              get_default_displacement_distance)

    supercell_matrix = np.zeros((3, 3))
    np.fill_diagonal(supercell_matrix, dim)

    if mode == 'chimes':
        mode = 'dftbp'

    # reading structure file
    cell_info = collect_cell_info(supercell_matrix,
                                  interface_mode=mode,
                                  return_dict=True)
    units = get_default_physical_units(interface_mode=mode)

    # Initialize phonopy
    phonon = Phonopy(cell_info['unitcell'],
	             cell_info['supercell_matrix'],
                     primitive_matrix=cell_info['primitive_matrix'],
                     factor=units['factor'],
                     calculator=cell_info['interface_mode'])

    # Create constant amplitude displacements
    d = get_default_displacement_distance(phonon.calculator)
    phonon.generate_displacements(distance=d)
    supercell = phonon.supercell
    supercells = phonon.supercells_with_displacements
    info = cell_info['optional_structure_info']
    additional_info = {'supercell_matrix': phonon.supercell_matrix}
    write_supercells_with_displacements(interface_mode=phonon.calculator,
                                        supercell=supercell,
                                        cells_with_disps=supercells,
                                        optional_structure_info=info,
                                        additional_info=additional_info)
    
    phonon.save()
    
def organize_folders(mode):
    """Finds supercell displacement files and formats the name. 
        Creates directory with supercell displacement number, 
    and moves supercell displacement file into created directory.

    Args:
        mode (str): Calculator used for task. Options are 'dftbp', 'chimes', 'vasp', or 'castep'. 
    """

    if mode == 'dftbp':
        for filename in os.listdir('.'):
            if filename.startswith('geo.genS-'):
                dir = filename[9:12]
                mkdir(dir)
                move(filename, '{}/geo_end.gen' .format(dir))

    if mode == 'chimes':
        for filename in os.listdir('.'):
            if filename.startswith('geo.genS-'):
                dir = filename[9:12]
                mkdir(dir)
                move(filename, '{}/geo_end.gen' .format(dir))
 
    if mode == 'vasp':
        for filename in os.listdir('.'):
            if filename.startswith('POSCAR-'):
                dir = filename[7:10]
                mkdir(dir)
                move(filename, '{}/POSCAR' .format(dir))

    if mode == 'castep':
        for filename in os.listdir('.'):
            if filename.startswith('supercell-'):
                dir = filename[10:13]
                mkdir(dir)
                move(filename, '{}/supercell.cell' .format(dir))


def calculate_forces(kforce, mode, dir):
    """Runs single point energy calculation.

    Args:
        kforce (list): Number of k points for force calculations.
        mode (str): Calculator used for task. Options are 'dftbp', 'chimes', 'vasp', or 'castep'. 
        dir (str): Directory to change to and run calculator in.
    """
    with chdir(dir):
        if isdone('forces'):
            return
        else:
            if mode == 'chimes':
                from ase.calculators.dftb import Dftb
                from dcs.chimes import run_md_input
                run_md_input('../../')
                calculator = Dftb(kpts=kforce,
                                  Hamiltonian_ChIMES='Yes',
                                  Hamiltonian_SCC='Yes',
                                  Hamiltonian_SCCTolerance=1e-7,
                                  Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}' .format(T=5),
                                  Hamiltonian_MaxAngularMomentum_='',
                                  Hamiltonian_MaxAngularMomentum_C='p',
                                  Hamiltonian_MaxAngularMomentum_O='p',
                                  Hamiltonian_MaxAngularMomentum_H='s',
                                  Hamiltonian_MaxAngularMomentum_N='p',
                                  Hamiltonian_MaxAngularMomentum_S='d',
                                  Analysis_='',
                                  Analysis_CalculateForces='Yes',
                                  Options_WriteResultsTag='Yes')
                calculator.write_dftb_in(filename='dftb_in.hsd')
                os.system('dftb+ 1>> forces.out 2>> forces.err')

            if mode == 'dftbp':
                from ase.calculators.dftb import Dftb
                calculator = Dftb(kpts=kforce,
                                  Hamiltonian_SCC='Yes',
                                  Hamiltonian_SCCTolerance=1e-7,
                                  Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}' .format(T=5),
                                  Hamiltonian_MaxAngularMomentum_='',
                                  Hamiltonian_MaxAngularMomentum_C='p',
                                  Hamiltonian_MaxAngularMomentum_O='p',
                                  Hamiltonian_MaxAngularMomentum_H='s',
                                  Hamiltonian_MaxAngularMomentum_N='p',
                                  Hamiltonian_MaxAngularMomentum_S='d',
                                  Analysis_='',
                                  Analysis_CalculateForces='Yes')
                calculator.write_dftb_in(filename='dftb_in.hsd')
                os.system('dftb+ 1>> forces.out 2>> forces.err')
                
            if mode == 'vasp':    
                from ase.calculators.vasp import Vasp
                atoms = read('POSCAR')
                calculator = Vasp(kpts=kforce,
                                  prec='Accurate',
                                  encut=520,
                                  ibrion=-1,
                                  ediff=1e-8,
                                  ismear=0,
                                  sigma=0.05,
                                  nwrite=1,
                                  ncore=16,
                                  lreal=False,
                                  lcharg=False,
                                  lwave=False,
                                  xc='pbe',
                                  gamma=True)
                calculator.calculate(atoms)
                
            if mode == 'castep':
                import ase.calculators.castep

                atoms = read('supercell.cell')
                calculator = ase.calculators.castep.Castep(kpts={'size':kforce, 'gamma':True})
                directory = '../' + dir
                calculator._export_settings = True
                calculator._directory = directory
                calculator._rename_existing_dir = False
                calculator._export_settings = True
                calculator._label = 'phonons'
                calculator.param.task = 'SinglePoint'
                calculator.param.xc_functional = 'PBE'
                calculator.param.cut_off_energy = 800
                calculator.param.elec_energy_tol = 1e-8
                calculator.param.num_dump_cycles = 0
                calculator.param.devel_code = 'PARALLEL: bands=4 kpoints=1 gvectors=64 :ENDPARALLEL'
                calculator.param.opt_strategy = 'speed'

                calculator.calculate(atoms)
                
            done('forces')
                
def multi_forces(kforce, mode, mpi=False):
    """Calls calculate_forces function using parallel processing,
        calculates forces for specified mode. 

    Args:
        kforce (list): Number of k points for force calculations.
        mode (str): Calculator used for task. Options are 'dftbp', 'chimes', 'vasp', or 'castep'. 
        mpi (bool, optional): Not currently implemented. Defaults to False. 
    """
    from functools import partial
    command = partial(calculate_forces, kforce, mode)
    
    dirlist = np.array(sorted([x.name for x in os.scandir() if x.is_dir()]))

    if mode == 'vasp' or 'castep':
        for dir in dirlist:
            command(dir)
    else:
        if mpi:
            from mpi4py.futures import MPIPoolExecutor
            with MPIPoolExecutor(max_workers=64, main=False) as executor:
                executor.map(command, dirlist)
        else:
            from multiprocessing import Pool        
            with Pool(processes=64) as pool:
                pool.map(command, dirlist)

    
def calculate_mesh(mesh, mode):
    """Creates Phonopy force sets for specified calculator and runs mesh sampling phonon calculation. 

    Args:
        mesh (list): Uniform meshes for each axis. 
        mode (str): Calculator used for task. Options are 'dftbp', 'chimes', 'vasp', or 'castep'.
    """
    import phonopy
    from phonopy.cui.create_force_sets import create_FORCE_SETS
    
    dirlist = sorted([x.name for x in os.scandir() if x.is_dir()])
    
    if mode == 'dftbp':
        filenames = [x + '/results.tag' for x in dirlist]
    if mode == 'chimes':
        mode = 'dftbp'
        filenames = [x + '/results.tag' for x in dirlist]
    if mode == 'vasp':
        filenames = [x + '/vasprun.xml' for x in dirlist]
    if mode == 'castep':
        filenames = [x + '/phonons.castep' for x in dirlist]

    create_FORCE_SETS(interface_mode=mode,
                      force_filenames=filenames,
                      disp_filename='phonopy_params.yaml')
    
    phonon = phonopy.load("phonopy_params.yaml")
    phonon.run_mesh(mesh, with_eigenvectors=True)
    phonon.write_yaml_mesh()


def phonons(dim=[4, 4, 4], kforce=[1, 1, 1], mesh=[8, 8, 8], calc='dftbp'):
    """Runs phonon supercell displacement calculations, populates 2-phonons folder with results.

    Args:
        dim (list, optional): Dimensions of the supercell. Defaults to [4, 4, 4].
        kforce (list, optional): Number of k points for force calculations. Defaults to [1, 1, 1].
        mesh (list, optional): Uniform meshes for each axis. Defaults to [8, 8, 8].
        calc (str, optional): Calculator used for task. Options are 'dftbp', 'chimes', 'vasp', or 'castep'. Defaults to 'dftbp'.

    Raises:
        NotImplementedError: Raised if 'calc' specified is not available. 
    """
    folder = os.getcwd()
    mkdir(folder + '/2-phonons')
    if calc == 'dftbp':
        copyfile(folder + '/1-optimization/geo_end.gen', folder + '/2-phonons/geo.gen')
    elif calc == 'vasp':
        copyfile(folder + '/1-optimization/CONTCAR', folder + '/2-phonons/POSCAR')
    elif calc == 'chimes':
        copyfile(folder + '/1-optimization/geo_end.gen', folder + '/2-phonons/geo.gen')
    elif calc == 'castep':
        unitcell = read(folder + '/1-optimization/relax.geom')
        write(folder + '/2-phonons/unitcell.cell', unitcell, positions_frac=True)

    else:
        raise NotImplementedError('{} calculator not implemented' .format(calc))
    
    with chdir(folder + '/2-phonons'):
        with out('phonons'):
            if isdone('phonons'):
                return
            else:
                generate_supercell(dim, calc)
                organize_folders(calc)
                multi_forces(kforce, calc)
                calculate_mesh(mesh, calc)
                done('phonons')

if __name__ == '__main__':
    phonons()