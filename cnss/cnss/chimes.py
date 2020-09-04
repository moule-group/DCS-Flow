import os
import numpy as np
from shutil import copyfile

default_atmmasses = {'C':12.0107,
                     'H':1.00784}

class CLICommand:
    'Chebyshev Interaction Model for Efficient Simulation (ChIMES)'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('--b2',
            help='order of 2 body interactions (integer)',
            default=12)
        add('--b3',
            help='order of 3 body interactions (integer)',
            default=8)
        add('--temp',
            help='Set the temperature in Kelvin',
            default=5)

    @staticmethod
    def run(args):
        chimes(args.b2, args.b3, args.temp)
        

def chimes_done():
    if os.path.exists('params.txt'):
        with open('params.txt') as f:
            if any('ENDFILE' in line for line in f.readlines()):
                return True
            else:
                return False
    else:
        return False
        
def dftb_fmatch_input(T):
    from ase.io import Trajectory
    from ase.calculators.dftb import Dftb
    from ase.units import Hartree, Bohr, GPa, mol, kcal

    traj = Trajectory('dft.traj')

    for frame in traj:
        n = frame.get_global_number_of_atoms()
        cell = frame.get_cell().lengths()
        cell = ' ' .join(map(str, cell))
        symbols = frame.get_chemical_symbols()
        positions = frame.get_positions()

        dft_forces = frame.get_forces()
        dft_stress = frame.get_stress()
        dft_energy = frame.get_total_energy()

        calc = Dftb(kpts=(1,1,1),
                    Hamiltonian_SCC='Yes',
                    Hamiltonian_MaxAngularMomentum_='',
                    Hamiltonian_MaxAngularMomentum_C='p',
                    Hamiltonian_MaxAngularMomentum_H='s',
                    Hamiltonian_Filling='Fermi {{Temperature [K] = {T} }}' .format(T=T))
        calc.calculate(frame, properties=['energy', 'forces', 'stress'])

        dftb_forces = calc.results['forces']
        dftb_stress = calc.results['stress']
        dftb_energy = calc.results['energy']

        diff_forces = (dft_forces - dftb_forces) / Hartree * Bohr

        diff_stress = -(dft_stress - dftb_stress)[:3] / GPa
        diff_stress = ' ' .join(map(str, diff_stress))

        diff_energy = (dft_energy - dftb_energy) * mol / kcal

        with open('dft-dftb.xyzf', 'a') as file:
            file.write('{} \n' .format(n))
            file.write('{} {} {} \n' .format(cell, diff_stress, diff_energy[0]))
            for s,p,f in zip(symbols, positions, diff_forces):
                p = ' ' .join(map(str, p))
                f = ' ' .join(map(str, f))
                file.write('{} {} {} \n' .format(s,p,f))

    nframes = len(traj)
    setsymbols = set(symbols)
    smax = round(min(frame.get_cell().lengths()) / 2, 2)
    return nframes, setsymbols, smax 

def rdf(smax, pair):
    from ase.io import Trajectory
    from ase.geometry.analysis import Analysis

    traj = Trajectory('dft.traj')
    ana = Analysis(traj[-1])
    rdf = ana.get_rdf(smax - 0.5, 100, elements=pair, return_dists=True)
    g = rdf[0][0]
    r = rdf[0][1]

    # choose rmin as the first position where g is greater than zero
    i = np.where(g > 0)
    rmin = r[i[0][0]]

    # smooth g function
    box_pts = 15
    box = np.ones(box_pts)/box_pts
    smoothg = np.convolve(g, box, mode='same')

    # choose rmax as first local minimum of smooth g
    ilocal = np.r_[True, smoothg[1:] < smoothg[:-1]] & np.r_[smoothg[:-1] < smoothg[1:], True]
    imax = (np.where(ilocal == True))[0][0]
    rmax = r[imax]

    # choose max value of g as morse lambda factor
    mlambda = max(g)

    # round up or down
    mlambda = round(mlambda, 2)
    rmin = float(str(rmin)[:3])
    rmax = float(str(rmax)[:3])
    
    return mlambda, rmin, rmax

def fm_setup_input(nframes, b2, b3, setsymbols, smax):
    if os.path.exists('../fm_setup.in'):
        copyfile('../fm_setup.in', 'fm_setup.in')
    else:
        with open('fm_setup.in', 'w') as f:
            f.write('####### CONTROL VARIABLES ####### \n'
                    '\n'
                    '# TRJFILE # \n'
                    '        dft-dftb.xyzf \n'
                    '# WRAPTRJ # \n'
                    '        false \n'
                    '# NFRAMES # \n'
                    '        {} \n'
                    '# NLAYERS # \n'
                    '        1 \n'
                    '# FITCOUL # \n'
                    '        false \n'
                    '# FITSTRS # \n'
                    '        false \n'
                    '# FITENER # \n'
                    '        false \n'
                    '# FITPOVR # \n'
                    '        false \n'
                    '# PAIRTYP # \n'
                    '        CHEBYSHEV {} {} 0 -1 1 \n'
                    '# CHBTYPE # \n'
                    '        MORSE \n'
	            '\n' .format(nframes, b2, b3))
            f.write('####### TOPOLOGY VARIABLES ####### \n'
                    '\n'
                    '# NATMTYP # \n'
                    '        {} \n'
                    '\n' .format(len(setsymbols)))
            f.write('# TYPEIDX #     # ATM_TYP #     # ATMCHRG #     # ATMMASS # \n')
            for idx, s in enumerate(setsymbols):
                f.write('{}               {}                0              {} \n'
                        .format(idx+1, s, default_atmmasses[s]))
            f.write('\n')
            f.write('# PAIRIDX #     # ATM_TY1 #     # ATM_TY1 #     # S_MINIM #     # S_MAXIM #     # S_DELTA #     # MORSE_LAMBDA #        # USEOVRP #     # NIJBINS #     # NIKBINS #     # NJKBINS # \n')

            from itertools import combinations_with_replacement
            pairs = list(combinations_with_replacement(setsymbols, 2))
            for idx, pair in enumerate(pairs):
                mlambda, rmin, rmax = rdf(smax, pair)
                f.write('{}               {}               {}               {}               {}             0.01              {}                    false           0               0               0 \n'
                        .format(idx+1, pair[0], pair[1], rmin, rmax, mlambda))
            f.write('\n'
                    '# FCUTTYP # \n'
                    '        TERSOFF 0.5 \n'
                    '\n'
                    '# ENDFILE #')

def run_md_input():
    if os.path.exists('../../run_md.in'):
        copyfile('../../run_md.in', 'run_md.in')
    elif os.path.exists('run_md.in'):
        return
    else:
        with open('run_md.in', 'w') as f:
            f.write('################################### \n'
                    '#### GENERAL CONTROL VARIABLES #### \n'
                    '################################### \n'
                    '\n'
                    '\n'
                    '# RNDSEED # ! Seed. If not specified, default value 123457 is used \n'
	            '        12357 \n'
                    '# TEMPERA # ! In K \n'
	            '        5.0 \n'
                    '# CMPRFRC # ! Compare computed forces against a set of input forces? ...If true, provide name of the file containing the forces for comparison \n'
	            '        false \n'
                    '# TIMESTP # ! In fs \n'
	            '        0.5 \n'
                    '# N_MDSTP # ! Total number of MD steps \n'
	            '        10000 \n'
                    '# NLAYERS # ! x,y, and z supercells.. small unit cell should have >= 1 \n'
	            '        1 \n'
                    '# USENEIG # ! Use a neighbor list? HIGHLY reccommended when NLAYERS > 0 \n'
	            '        true \n'
                    '# PRMFILE # ! Parameter file (i.e. params.txt) \n'
	            '        params.txt \n'
                    '# CRDFILE # ! Coordinate file (.xyz) or force file (.xyzf) \n'
	            '        dftb.xyz \n'
                    '# TRAJEXT # ! coordinate file type \n'
	            '        GEN \n'
	            '\n'
	            '\n'
                    '################################### \n'
                    '####    SIMULATION  OPTIONS    #### \n'
                    '################################### \n'
                    '\n'
                    '# VELINIT # (options are READ or GEN) \n'
	            '        GEN \n'
                    '# CONSRNT # (options are HOOVER <hoover time> or VELSCALE <scale freq> \n'
	            '        NVT-MTK HOOVER 10 \n'
                    '# PRSCALC # (options are ANALYTICAL or NUMERICAL) \n'
	            '        ANALYTICAL \n'
	            '\n'
                    '################################### \n'
                    '####      OUTPUT  CONTROL      #### \n'
                    '################################### \n'
	            '\n'
                    '# WRPCRDS # \n'
	            '        false \n'
                    '# FRQDFTB # ! Frequency to output the DFTB gen file \n'
	            '        10 \n'
                    '# FRQENER # ! Frequency to output energies \n'
	            '        1 \n'
                    '# PRNTFRC # ! Print computed forces? Forces are printed to force_out.txt \n'
	            '        true 20 \n'
	            '\n'
                    '# ENDFILE #')
    
    
def lsq():
    os.system('chimes_lsq fm_setup.in')
    os.system('lsq2.py --algorithm=lassolars > params.txt')
    
    
def chimes(b2=12, b3=8, T=5):
    from cnss import mkdir, chdir, out
    import glob

    folder = os.getcwd()
    trajfile = glob.glob(folder + '/1_1-molecular_dynamics/*.traj')[0]
    mkdir(folder + '/1_2-chimes')
    copyfile(trajfile, folder + '/1_2-chimes/dft.traj')
    with chdir(folder + '/1_2-chimes'):
        with out('chimes'):
            if chimes_done():
                return
            else:
                nframes, setsymbols, smax = dftb_fmatch_input(T)
                fm_setup_input(nframes, b2, b3, setsymbols, smax)
                lsq()
    

if __name__ == '__main__':
    chimes()
