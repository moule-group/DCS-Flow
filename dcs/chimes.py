import os
import numpy as np
from shutil import copyfile

default_atmmasses = {'C':12.0107,
                     'H':1.00784,
                     'O':15.999,
                     'N':14.0067,
                     'S':32.065}

class CLICommand:
    'Chebyshev Interaction Model for Efficient Simulation (ChIMES)'

    @staticmethod
    def add_arguments(parser):
        """Sets up command line to run Chimes i.e. recognize arguments and commands.

        Args:
            parser (argparse): Arguments to be added. 
        """
        add = parser.add_argument
        add('--trajfile',
            help='path to DFT-MD trajectory file',
            default=None)
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
        """Runs chimes function using command line arguments. 

        Args:
            args (argparse): Command line arguments added to parser using the function add_arguments.
        """
        chimes(args.trajfile, args.b2, args.b3, args.temp)
        
        
def dftb_fmatch_input(T, frame):
    """Runs DFTB force calculations using frame from DFT;
        calculates the force difference between DFT and DFTB.
        Creates a list with the chemical symbols, atomic positions, and force differences. 

    Args:
        T (int): Temperature for DFTB simulation.
        frame (list): ASE atoms object with atomic positions, forces, energies, etc from DFT calculation.

    Returns:
        list: Contains the chemical symbols, atom positions, and difference in DFTB and DFT forces.
    """
    from ase.calculators.dftb import Dftb
    from ase.units import Hartree, Bohr, GPa, mol, kcal
    from dcs import mktempdir, chdir

    with mktempdir() as temp_dir:
        with chdir(temp_dir):
            n = frame.get_global_number_of_atoms()
            cell = frame.get_cell().lengths()
            cell = ' ' .join(map(str, cell))
            symbols = frame.get_chemical_symbols()
            positions = frame.get_positions()
    
            dft_forces = frame.get_forces()
            # dft_stress = frame.get_stress()
            dft_energy = frame.get_total_energy()

            calc = Dftb(label='dftb',
                        kpts=(1,1,1),
                        Hamiltonian_SCC='Yes',
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}' .format(T=T))
            # calc.calculate(frame, properties=['energy', 'forces', 'stress'])
            calc.calculate(frame, properties=['energy', 'forces'])

            # saving output from DFTB+ to chimes.out
            with open('dftb.out') as f:
                print(f.read())

            dftb_forces = calc.results['forces']
            # dftb_stress = calc.results['stress']
            dftb_energy = calc.results['energy']

            diff_forces = (dft_forces - dftb_forces) / Hartree * Bohr

            # diff_stress = -(dft_stress - dftb_stress)[:3] / GPa
            # diff_stress = ' ' .join(map(str, diff_stress))

            diff_energy = (dft_energy - dftb_energy) * mol / kcal

            frame_fmatch = []

            frame_fmatch.append('{} \n' .format(n))
            # frame_fmatch.append('{} {} {} \n' .format(cell, diff_stress, diff_energy))
            frame_fmatch.append('{} {} \n' .format(cell, diff_energy))
            for s,p,f in zip(symbols, positions, diff_forces):
                p = ' ' .join(map(str, p))
                f = ' ' .join(map(str, f))
                frame_fmatch.append('{} {} {} \n' .format(s,p,f))

    return frame_fmatch

def multi_fmatch(T):
    """Runs dftb_fmatch_input command for each trajectory, n processes at a time.
        Returns tuple with trajectory information. 

    Args:
        T (int): Temperature for DFTB simulation. 

    Returns:
        tuple: (nframes, setsymbols, smax)
                nframes (int) - Number of ASE atoms object.
                setsymbols (set) - Chemical symbols. 
                smax (float) - Half of the minimum cell length.
    """
    from ase.io import Trajectory
    from multiprocessing import Pool
    from functools import partial

    traj = Trajectory('dft.traj')

    command = partial(dftb_fmatch_input, T)
    with Pool(processes=4) as pool:
        results = pool.map(command, list(traj))

    with open('dft-dftb.xyzf', 'a') as file:
        for result in results:
            for line in result:
                file.write(line)

    nframes = len(traj)
    setsymbols = set(traj[0].get_chemical_symbols())
    smax = round(min(traj[0].get_cell().lengths()) / 2, 2)

    return nframes, setsymbols, smax 

def rdf(smax, pair):
    """Finds the radial distribution function for elements defined in pair. 
       Analyses the RDF, and determines the min and max distances in a atomic pair interaction.
       Returns a tuple with mlambda, rmin and rmax.

    Args:
        smax (float): Half of the minimum cell length.
        pair (list): List with two chemical symbols of a pair interaction, e.g., ['C', 'H']

    Returns:
        tuple: (mlambda, rmin, rmax)
            mlambda (float) - Morse lambda factor
            rmin (float) - Minimum distance considered in a atomic pair interaction.
            rmax (float) - Maximum distance considered in a atomic pair interaction.
    """
    from ase.io import Trajectory
    from ase.geometry.analysis import Analysis

    traj = Trajectory('dft.traj')
    ana = Analysis(list(traj))

    rdf = ana.get_rdf(smax-0.4, 100, elements=pair, return_dists=True)
    r = rdf[0][1]

    # get average g function                                                                            
    g = []
    for i in range(len(traj)):
        g.append(rdf[i][0])
    aveg = np.average(g, axis=0)

    import matplotlib.pyplot as plt
    plt.figure(pair[0] + '-' + pair[1])
    plt.plot(r, aveg)
    plt.savefig(pair[0] + '-' + pair[1] + '.png', dpi=300, bbox_inches='tight', pad_inches=0)
    
    # choose rmin as the first position where aveg is greater than zero          
    i = np.where(aveg > 0)
    rmin = r[i[0][0]]

    # smooth g function
    box_pts = 5
    box = np.ones(box_pts)/box_pts
    smoothg = np.convolve(aveg, box, mode='same')
    
    # choose rmax as first zero after peak or local minimum of smooth g
    try:
        ilocal = np.where(np.r_[True, smoothg[1:] < smoothg[:-1]] & np.r_[True, smoothg[1:] == 0])
        imax = ilocal[0][1]
        rmax = r[imax]
    except:
        try:
            ilocal = np.r_[True, smoothg[1:] < smoothg[:-1]] & np.r_[smoothg[:-1] < smoothg[1:], True]
            imax = (np.where(ilocal == True))[0][0]
            rmax = r[imax]
        except:
            if smax < 2 * rmin:
                rmax = smax - 0.4
            else:
                rmax = 2 * rmin

    # choose position of max value of aveg as morse lambda factor                                       
    igmax = np.argmax(aveg)
    mlambda = r[igmax]

    # round up or down                                                                                 
    mlambda = round(mlambda, 2)
    rmin = float(str(rmin)[:3])
    rmax = float(str(rmax)[:3])

    return mlambda, rmin, rmax


def fm_setup_input(nframes, b2, b3, setsymbols, smax):
    """Writes force match setup input file (fm_setup.in) with given arguments.

    Args:
        nframes (int): Number of trajectories.
        b2 (int): Second body order of Chebyshev polynomial. 
        b3 (int): Third body order of Chebyshev polynomial.
        setsymbols (set): Chemical symbols. 
        smax (float): Half of the minimum cell length.
    """
    if os.path.exists('../../fm_setup.in'):
        copyfile('../../fm_setup.in', 'fm_setup.in')
    elif os.path.exists('../fm_setup.in'):
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

def run_md_input(folder):
    """Writes MD input file (run_md.in) if it does not exist in folder.

    Args:
        folder (str): Directory that contains params.txt file.

    Raises:
        Exception: Raises exception if there is no params.txt file in folder.
    """
    if os.path.exists(folder + '/params.txt'):
        copyfile(folder + '/params.txt', 'params.txt')
    else:
        raise Exception('you need a ChIMES parameters file called params.txt')
    
    if os.path.exists(folder + '/run_md.in'):
        copyfile(folder + '/run_md.in', 'run_md.in')
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
	            '        1.0 \n'
                    '# N_MDSTP # ! Total number of MD steps \n'
	            '        10000 \n'
                    '# NLAYERS # ! x,y, and z supercells.. small unit cell should have >= 1 \n'
	            '        0 \n'
                    '# USENEIG # ! Use a neighbor list? HIGHLY reccommended when NLAYERS > 0 \n'
	            '        false \n'
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
    """Runs chimes_lsq binary and creates A and B matrices (A.txt, B.txt).
       Runs lsq fitting process and writes ChIMES coefficients to params.txt.
    """
    os.system('chimes_lsq fm_setup.in >> chimes.out')
    os.system('lsq2.py --algorithm=lassolars > params.txt')
    
    
def chimes(trajfile=None, b2=12, b3=8, T=5):
    """Calculates force difference between DFT and DFTB (training set), 
       and fits the Chebyshev polynomials coefficients.
       Creates 3-chimes folder (inside 0-train directory) and writes params.txt file.

    Args:
        trajfile (list, optional): Trajectory file (ist of atoms objects) output from md simulation. Defaults to None.
        b2 (int, optional): Second body order of Chebyshev polynomial. Defaults to 12.
        b3 (int, optional): Third body order of Chebyshev polynomial. Defaults to 8.
        T (int, optional): Temperature for simulation in Kevin. Defaults to 5.
    """
    from dcs import mkdir, chdir, out, done, isdone
    import glob

    folder = os.getcwd()
    if trajfile is None:
        if os.path.isdir(folder + '/2-molecular_dynamics'):
            trajfile = glob.glob(folder + '/2-molecular_dynamics/*.traj')[0]
        else:
            trajfile = glob.glob('*.traj')[0]
            
    mkdir(folder + '/3-chimes')
    copyfile(trajfile, folder + '/3-chimes/dft.traj')
    with chdir(folder + '/3-chimes'):
        with out('chimes'):
            if isdone('chimes'):
                return
            else:
                nframes, setsymbols, smax = multi_fmatch(T)
                fm_setup_input(nframes, b2, b3, setsymbols, smax)
                lsq()
                done('chimes')

if __name__ == '__main__':
    chimes()
