import os
import glob
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.io import read, write
from dcs import mkdir, chdir, out, done, isdone

class CLICommand:
    'Molecular dynamics with constant temperature'

    @staticmethod
    def add_arguments(parser):
        """Sets up command line to run molecular dynamics with constant temperature i.e. recognize arguments and commands.

        Args:
            parser (argparse): Arguments to be added. 
        """
        add = parser.add_argument
        add('--optgeo',
            help='path to DFT optimized geometry file',
            default=None)
        add('--calc',
            help='Calculator used for moculecular dynamics simulation. Options are vasp and castep',
            default='vasp')
        add('--temp',
            help='Set the temperature in Kelvin',
            default=300,
            type=int)
        add('--md_size',
            help='Size of supercell for md, e. g., 1 1 1',
            default=[1, 1, 1],
            nargs=3,
            type=int)
        add('--steps',
            help='Set the number of steps for MD simulation',
            default=5000,
            type=int)
        add('--time_step',
            help='Set the time step in fs',
            default=1,
            type=int)
        add('--dump_interval',
            help='Set the interval to save the frame on the trajectory file',
            default=100,
            type=int)


    @staticmethod
    def run(args):
        """Runs md function using command line arguments. 

        Args:
            args (argparse): Command line arguments added to parser using the function add_arguments.
        """
        md(args.optgeo, args.calc, args.temp, args.md_size,
           args.steps, args.time_step, args.dump_interval)


def run_vasp_md(atoms, T, steps, time_step, dump_interval):
    """Runs vasp md calculation on atoms for specified intervals.

    Args:
        atoms (list): Atoms object involved in calc from ASE.
        T (int): Initial and final simulation temperature.
        steps (int): Maximum number of ionic steps, defines the total simulation time.
        time_step (int): md time step in fs.
        dump_interval (int): Step size of frames to be saved in the trajectory file.
    """
    if isdone('md'):
        return
    else:
        from ase.calculators.vasp import Vasp
        atoms.calc = Vasp(encut=520,
                          prec='Normal',
                          algo='Fast', # electronic minimisation algotithm
                          lreal='Auto', # operators in real space
                          ismear=0, # Gaussian smearing
                          sigma=0.05, # width of smearing in eV
                          # ismear=-1, # Fermi smearing
                          # sigma=T*units.kB, # width of smearing in eV
                          isym=0, # no symmetry usage                                                   
                          ibrion=0, # MD run
                          potim=time_step, # md time step                                              
                          nsw=steps, # maximum number of ionic steps                                    
                          tebeg=T, # initial temperature                                                
                          teend=T, # final temperature                                                  
                          mdalgo=2, # Nose-Hoover thermostat                                            
                          isif=2, # NVT ensemble                                                        
                          smass=0, # Canonical (Nose-Hoover) thermostat                                 
                          ediff=1e-6, # global break condition for the electronic SC-loop               
                          nwrite=1, # how much will be written to the OUTCAR file                      
                          ncore=16, # number of bands that are treated in parallel
                          lcharg=False, # charge densities are not written                              
                          lwave=False, # wavefunctions are not written                                  
                          xc='pbe')
        atoms.get_potential_energy()
        vasptraj = read('OUTCAR', index=slice(0, steps, dump_interval))
        write('md.traj', vasptraj)
        done('md')

def run_castep_md(atoms, T, steps, time_step, dump_interval):
    """Runs castep md calculation on atoms for specified intervals.

    Args:
        atoms (list): Atoms object involved in calc from ASE.
        T (int): Initial and final simulation temperature.
        steps (int): Maximum number of ionic steps, defines the total simulation time.
        time_step (int): md time step in fs.
        dump_interval (int): Step size of frames to be saved in the trajectory file.
    """
    if isdone('md'):
        return
    else:
        import ase.calculators.castep
        calculator = ase.calculators.castep.Castep(kpts={'size':[1,1,1], 'gamma':True})
        directory = '../2-molecular_dynamics'
        calculator._export_settings = True
        calculator._directory = directory
        calculator._rename_existing_dir = False
        calculator._export_settings = True
        calculator._label = 'md'
        calculator._set_atoms = True
        calculator.param.task = 'MolecularDynamics'
        calculator.param.xc_functional = 'PBE'
        calculator.param.cut_off_energy = 800
        calculator.param.opt_strategy = 'speed'
        calculator.param.calculate_stress = True
        calculator.param.md_sample_iter = dump_interval
        calculator.param.popn_calculate = False
        calculator.param.num_dump_cycles = 0
        calculator.param.md_num_iter = steps
        calculator.param.md_delta_t = str(time_step) + ' fs'
        calculator.param.md_ensemble = 'NVT'
        calculator.param.md_temperature = str(T) + ' K'
        calculator.param.md_thermostat = 'nose-hoover'
        calculator.param.md_elec_energy_tol = 1e-6
        calculator.param.devel_code = 'PARALLEL: bands=4 kpoints=1 gvectors=64 :ENDPARALLEL'
        calculator.cell.fix_com = True
        calculator.cell.fix_all_cell = True
        calculator.cell.symmetry_generate = True
        atoms.set_calculator(calculator)
        atoms.get_potential_energy()
        casteptraj = read('md.md', index=slice(0, steps, dump_interval))
        write('md.traj', casteptraj)
        done('md')


def md(optgeo=None, calc='vasp', T=300, md_size=[1,1,1],
       steps=5000, time_step=1, dump_interval=100):
       """Runs molecular dynamics simulation using vasp or castep
       (Creates 2-molecular_dynamics folder inside 0-train)

       Args:
       optgeo (NoneType, optional): Optimized geometry file, only true if optgeo defined. 
       Defaults to None.
       calc (str, optional): Specifies calculator. Options are 'vasp' or 'castep'. Defaults to 'vasp'.
       T (int, optional): Simulation temperature. Defaults to 300.
       md_size (list, optional): Size of supercell. Defaults to [1,1,1].
       steps (int, optional): Maximum number of ionic steps. Defaults to 5000.
       time_step (int, optional): Md time step in fs. Defaults to 1. 
       dump_interval (int, optional): Step size of frames to be saved in the trajectory file. 
       Defaults to 100. 

       Raises:
       NotImplementedError: If calculator other than vasp or castep specified.
       """

       folder = os.getcwd()

       if optgeo:
           atoms = read(optgeo)
       else:
           if calc == 'vasp':
               atoms = read(folder + '/1-optimization/CONTCAR')
           elif calc == 'castep':
               atoms = read(folder + '/1-optimization/relax.geom')
           else:
               raise NotImplementedError('{} calculator not implemented' .format(calc))
           

       atoms = atoms.repeat(md_size)
       mkdir(folder + '/2-molecular_dynamics')
       with chdir(folder + '/2-molecular_dynamics'):
           with out('md'):
               if calc=='vasp':
                   run_vasp_md(atoms, T, steps, time_step, dump_interval)
               elif calc=='castep':
                   run_castep_md(atoms, T, steps, time_step, dump_interval)
               else:
                   raise NotImplementedError('{} calculator not implemented' .format(calc))

if __name__ == '__main__':
    md()

