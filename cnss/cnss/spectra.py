import os
import numpy as np
from shutil import copyfile, move
from ase.io import read, write
from cnss import mkdir, chdir, out, done, isdone
from phonopy import Phonopy

################# STILL IN INITIAL STATE ########################

class CLICommand:
    'Simulate spectroscopy with supercell phonon method'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('--calc',
            help='Calculator used. Option is castep',
            default='castep')
        add('--type',
            help='Type of spectroscopy simulation. Options are ir and raman',
            default='ir')
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
        phonons(args.dim, args.kforce, args.mesh, args.calc, args.type)


def ir():
    import ase.calculators.castep
    
    atoms = read('supercell.cell')
    calculator = ase.calculators.castep.Castep()
    directory = '../' + dir
    calculator._export_settings = True
    calculator._directory = directory
    calculator._rename_existing_dir = False
    calculator._export_settings = True
    calculator._label = 'phonons'
    calculator.param.task = 'SinglePoint'
    calculator.param.xc_functional = 'PBE'
    calculator.param.cut_off_energy = 520
    calculator.param.elec_energy_tol = 1e-8
    calculator.param.num_dump_cycles = 0
    calculator.cell.kpoint_mp_grid = kforce
    
    calculator.calculate(atoms)

        
def raman():
    import ase.calculators.castep
    
    atoms = read('supercell.cell')
    calculator = ase.calculators.castep.Castep()
    directory = '../' + dir
    calculator._export_settings = True
    calculator._directory = directory
    calculator._rename_existing_dir = False
    calculator._export_settings = True
    calculator._label = 'phonons'
    calculator.param.task = 'SinglePoint'
    calculator.param.xc_functional = 'PBE'
    calculator.param.cut_off_energy = 520
    calculator.param.elec_energy_tol = 1e-8
    calculator.param.num_dump_cycles = 0
    calculator.cell.kpoint_mp_grid = kforce
    
    calculator.calculate(atoms)



def spectra(dim=[4, 4, 4], kforce=[1, 1, 1], mesh=[8, 8, 8], calc='castep', type='ir'):
    folder = os.getcwd()
    mkdir(folder + '/2-spectra')
    if calc == 'castep':
        unitcell = read(folder + '/1-optimization/relax.geom')
        write(folder + '/2-spectra/unitcell.cell', unitcell, positions_frac=True)

    else:
        raise NotImplementedError('{} calculator not implemented' .format(calc))
    
    with chdir(folder + '/2-spectra'):
        with out('spectra'):
            if isdone('spectra'):
                return
            else:
                if type == 'ir':
                    ir()
                elif type == 'raman':
                    raman()
                else:
                    raise NotImplementedError('{} not implemented' .format(calc))
                
                done('spectra')

if __name__ == '__main__':
    spectra()
