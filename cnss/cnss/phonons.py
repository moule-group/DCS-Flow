import os
import numpy as np
import glob
from shutil import copyfile, move
from ase.calculators.dftb import Dftb
from ase.io import read
from cnss import mkdir, chdir, out
import phonopy
from phonopy import Phonopy

class CLICommand:
    'Calculate phonons'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
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
        phonons(args.dim, args.kforce, args.mesh)


def generate_supercell(dim, d=0.01):
    from phonopy.interface.calculator import read_crystal_structure
    from phonopy.interface.calculator import write_supercells_with_displacements
    from phonopy.interface.calculator import get_default_physical_units

    mode='dftbp'

    unitcell, info = read_crystal_structure('geo.gen', interface_mode=mode)
    supercell_matrix = np.zeros((3, 3))
    np.fill_diagonal(supercell_matrix, dim)

    units = get_default_physical_units(interface_mode=mode)

    phonon = Phonopy(unitcell=unitcell,
                     supercell_matrix=supercell_matrix,
                     factor=units['factor'])
    phonon.generate_displacements(distance=d)
    supercell = phonon.supercell
    supercells = phonon.supercells_with_displacements
    
    write_supercells_with_displacements(interface_mode=mode,
                                        supercell=supercell,
                                        cells_with_disps=supercells,
                                        optional_structure_info=info,
                                        additional_info={})
    return phonon
    
def organize_folders():

    for filename in os.listdir('.'):
        if filename.startswith('geo.genS-'):
            dir = filename[9:12]
            mkdir(dir)
            move(filename, '{}/geo_end.gen' .format(dir))

            
def write_dftb_input(kforce):
    dirlist = sorted([x.name for x in os.scandir() if x.is_dir()])
    for dir in dirlist:
        with chdir(dir):
            atoms = read('{}' .format(glob.glob('*.gen')[0]))
            calc = Dftb(atoms=atoms,
                        kpts=kforce,
                        Hamiltonian_SCC='No',
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Analysis_='',
                        Analysis_CalculateForces='Yes',
                        Options_WriteResultsTag='Yes')
            atoms.set_calculator(calc)
            calc.write_dftb_in(filename='dftb_in.hsd')


def dftb_command(dir):
    with chdir(dir):
        if not os.path.exists('results.tag') or not os.path.getsize('results.tag') > 0:
            os.system('dftb+ 1>> forces.out 2>> forces.err')

        
def calculate_forces(mpi=False):
    
    dirlist = np.array(sorted([x.name for x in os.scandir() if x.is_dir()]))
    
    if mpi:
        from mpi4py.futures import MPIPoolExecutor
        with MPIPoolExecutor(max_workers=16, main=False) as executor:
            executor.map(dftb_command, dirlist)
    else:
        from multiprocessing import Pool        
        pool = Pool(processes=16)
        pool.map(dftb_command, dirlist)


def calculate_mesh(phonon, mesh):
    from phonopy.interface.calculator import get_force_sets

    dirlist = sorted([x.name for x in os.scandir() if x.is_dir()])
    
    mode = 'dftbp'
    if mode == 'dftbp':
        filenames = [x + '/results.tag' for x in dirlist]

    natom = phonon.supercell.get_number_of_atoms()
    forces = get_force_sets(interface_mode=mode,
                            num_atoms=natom,
                            num_displacements=len(filenames),
                            force_filenames=filenames)
    phonon.set_forces(forces)
    phonon.produce_force_constants()
    phonon.save()
    phonon.run_mesh(mesh, with_eigenvectors=True)
    phonon.write_yaml_mesh()


def phonons(dim=[4, 4, 4], kforce=[1, 1, 1], mesh=[8, 8, 8]):
    folder = os.getcwd()
    mkdir(folder + '/2-phonons')
    copyfile(folder + '/1-optimization/geo_end.gen', folder + '/2-phonons/geo.gen')
    
    with chdir(folder + '/2-phonons'):
        with out('phonons'):
            phonon = generate_supercell(dim)
            organize_folders()
            write_dftb_input(kforce)
            calculate_forces()
            calculate_mesh(phonon, mesh)

if __name__ == '__main__':
    phonons()
