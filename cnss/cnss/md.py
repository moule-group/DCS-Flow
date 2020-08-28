import os
import glob
from ase import units
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from cnss import mkdir, chdir, out

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
            default=300)
        add('--md_size',
            help='Size of supercell for md, e. g., 2 2 2',
            default=[2, 2, 2],
            nargs=3,
            type=int)


    @staticmethod
    def run(args):
        md(args.md_calc, args.temp, args.md_size)


def md_done(steps):
    if os.path.exists('md.out'):
        with open('md.out') as f:
            if len(f.readlines()) > steps:
                return True
            else:
                return False
    else:
        return False    
    
def print_energy(a):
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

    
def run_vasp_md(atoms, T):
    steps = 20
    if md_done(steps):
        return
    else:
        from ase.calculators.vasp import Vasp
        atoms.calc = Vasp(encut=520,
                          prec='Accurate',
                          nwrite=1,
                          npar=8,
                          lreal='Auto',
                          lcharg=False,
                          lwave=False,
                          xc='pbe')
        dyn = Langevin(atoms, 1 * units.fs, T * units.kB, 0.002)
        dyn.attach(print_energy, 1, atoms)
        traj = Trajectory('md.traj', 'w', atoms)
        dyn.attach(traj.write, interval=1)
        print_energy(atoms)
        dyn.run(steps)

    
def md(md_calc='vasp', T=300, md_size=[2,2,2]):
    folder = os.getcwd()
    trajfile = glob.glob(folder + '/1-optimization/*.traj')[0]
    
    mkdir(folder + '/1_1-molecular_dynamics')
    with chdir(folder + '/1_1-molecular_dynamics'):
        with out('md'):
            traj = Trajectory(trajfile)
            atoms = traj[-1]
            atoms = atoms.repeat(md_size)

            if md_calc=='vasp':
                run_vasp_md(atoms, T)
            else:
                raise NotImplementedError('{} calculator not implemented' .format(md_calc))

if __name__ == '__main__':
    md()
