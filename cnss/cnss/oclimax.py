import os
import argparse
from cnss import mkdir, chdir
from shutil import copyfile

class CLICommand:
    'Run oclimax simulation'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('--params',
            help='Parameters file, e. g., out.params',
            default=None)
        add('--task',
            help='0:inc approx. 1:coh+inc. 2:single-xtal Q-E. 3:single-xtal Q-Q',
            default=1)
        add('--e_unit',
            help='Energy unit [eu] (0:cm-1,1:meV,2:THz)',
            default=0)

    @staticmethod
    def run(args):
        oclimax(args.params, args.task, args.e_unit)


def write_params(task, e_unit):
    with open('out.params', 'w') as f:
        f.write('## General parameters \n'
                'TASK    =         {} # 0:inc approx. 1:coh+inc. 2:single-xtal Q-E. 3:single-xtal Q-Q\n'
                'INSTR   =         0  # 0:VISION 1:indirect traj 2:direct traj 3:Q-E or Q-Q mesh\n'
                'TEMP    =      0.00  # Temperature [K]\n'
                'E_UNIT  =         {} # Energy unit [eu] (0:cm-1,1:meV,2:THz)\n'
                'OUTPUT  =   0  # 0:standard, 1:restart, 2:SPE, 3:full, 4:DOS, 5:modes\n'
 
                '## Additional general parameters\n'
                'MAXO    =        10  # Maximum order of excitation\n'
                'CONV    =         2  # Start convolution from order=CONV (2 or 3)\n'
                'PHASE   =         0  # Phase factor of polarization vector (0 or 1)\n'
                'MASK    =         0  # Set 1 to apply mask on Q-E/Q-Q map (INSTR=3)\n'
                'ELASTIC =  -0.10E+01 -0.10E+01  # E Q, <0:no EL,0:cal res,>0:given res\n'
 
                '## E parameters\n'
                'MINE    =      8.00  # Energy range (minimum) to calculate [eu]\n'
                'MAXE    =   3500.00  # Energy range (maximum) to calculate [eu]\n'
                'dE      =      1.00  # Energy bin size [eu]\n'
                'ECUT    =     8.000  # Exclude modes below this cutoff energy [eu]\n'
                'ERES    =   0.25E+01  0.50E-02  0.10E-06  # E resolution coeff\n'
 
                '## Q parameters\n'
                'MINQ    =      0.50  # Q range (minimum) to calculate\n'
                'MAXQ    =     20.00  # Q range (maximum) to calculate\n'
                'dQ      =      0.50  # Q bin size\n'
                'QRES    =   0.10E+00 # Q resolution coeff (INSTR=3)\n'
 
                '## Instrument parameters\n'
                'THETA   =  135.0  45.0  # List of scattering angles [degree]\n'
                'Ef      =     32.00  # Final energy [eu] (INSTR=1)\n'
                'Ei      =   5000.00  # Incident energy [eu] (INSTR=2)\n'
                'L1      =     11.60  # L1 [m] for DGS (INSTR=2 or 3, ERES=0)\n'
                'L2      =      2.00  # L2 [m] for DGS (INSTR=2 or 3, ERES=0)\n'
                'L3      =      3.00  # L3 [m] for DGS (INSTR=2 or 3, ERES=0)\n'
                'dt_m    =      3.91  # dt_m [us] for DGS (INSTR=2 or 3, ERES=0)\n'
                'dt_ch   =      5.95  # dt_ch [us] for DGS (INSTR=2 or 3, ERES=0)\n'
                'dL3     =      3.50  # dL3 [cm] for DGS (INSTR=2 or 3, ERES=0)\n'
 
                '## Single crystal parameters\n'
                'HKL     =    0.0   0.0   0.0  # HKL (TASK=2 or 3)\n'
                'Q_vec   =    0.0   0.0   1.0  # Q vector dir (TASK=2 or 3)\n'
                'Q_vec_y =    1.0   0.0   0.0  # Q vector dir y-axis (TASK=3)\n'
                'MINQ_y  =      1.00  # Q range (minimum) y-axis (TASK=3)\n'
                'MAXQ_y  =      2.00  # Q range (maximum) y-axis (TASK=3)\n'
                'dQ_y    =      0.02  # Q bin size y-axis (TASK=3)\n'
                
                '## Wing parameters\n'
                'WING    =         0  # Wing calculation (0:no wing,1:isotropic,2:ST tensor)\n'
                'A_ISO   =    0.0350  # Isotropic A_external for wing calculation\n'
                'W_WIDTH =     150.0  # Energy width [eu] of initial wing)\n' .format(task, e_unit))
    

def run_oclimax(params):
    
    os.system('oclimax convert -yaml mesh.yaml -o > ocl.out')
    os.system('oclimax run out.oclimax {} >> ocl.out' .format(params))

def plot():
    import pandas as pd
    import matplotlib.pyplot as plt
    import glob

    file = glob.glob('*.csv')
    df = pd.read_csv(file[0], skiprows=4, header=None)
    E = df.iloc[:, 0]
    totback = df.iloc[:, 1]
    totfor = df.iloc[:, 2]
    int = (totback + totfor) / 2
    normint = int
    # normint = ((int - min(int[200:])) / (max(int[200:]) - min(int[200:])))

    plt.plot(E, normint)
    plt.xlabel('Energy (cm$^{-1}$)')
    plt.ylabel('Normalized intensity')
    # plt.xlim(0, 3500)
    # plt.ylim(0, 1)
    plt.savefig(file[0][:-4]+'.png', dpi=300, bbox_inches='tight', pad_inches=0)


def oclimax(params=None, task=1, e_unit=0):
    folder = os.getcwd()
    mkdir(folder + '/3-oclimax')
    copyfile(folder + '/2-phonons/mesh.yaml', folder + '/3-oclimax/mesh.yaml')

    with chdir(folder + '/3-oclimax'):
        if not params:
            write_params(task, e_unit)
            params = 'out.params'
        run_oclimax(params)
        plot()

if __name__ == '__main__':
    oclimax()
