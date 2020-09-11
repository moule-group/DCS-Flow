import os
import glob
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.io import read, write
from cnss import mkdir, chdir, out, done, isdone

class CLICommand:
    'Molecular dynamics with constant temperature'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('--md_calc',
            help='Calculator used for moculecular dynamics simulation. Options are vasp',
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


    @staticmethod
    def run(args):
        md(args.md_calc, args.temp, args.md_size)


def run_vasp_md(atoms, T):
    steps = 100
    time_step = 1 # in fs                                                                               
    dump_interval = 1
    if isdone('md'):
        return
    else:
        from ase.calculators.vasp import Vasp
        atoms.calc = Vasp(encut=550,
                          prec='Normal',
                          algo='Fast', # electronic minimisation algotithm                              
                          lreal='Auto', # operators in real space                                       
                          ismear=0, # gaussian smearing
                          sigma=0.1, # width of smearing in eV
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
                          npar=8, # number of bands that are treated in parallel                       
                          lcharg=False, # charge densities are not written                              
                          lwave=False, # wavefunctions are not written                                  
                          xc='pbe')
        atoms.get_potential_energy()
        vasptraj = read('OUTCAR', index=slice(0, steps, dump_interval))
        write('md.traj', vasptraj)
        done('md')
        

def md(md_calc='vasp', T=300, md_size=[1,1,1]):
    folder = os.getcwd()
    mkdir(folder + '/1_1-molecular_dynamics')
    with chdir(folder + '/1_1-molecular_dynamics'):
        with out('md'):
            atoms = read(folder + '/1-optimization/POSCAR')
            atoms = atoms.repeat(md_size)

            if md_calc=='vasp':
                run_vasp_md(atoms, T)
            else:
                raise NotImplementedError('{} calculator not implemented' .format(md_calc))

if __name__ == '__main__':
    md()


