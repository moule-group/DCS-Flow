# Davis Computational Spectroscopy workflow (DCS-Flow)

DCS-Flow is a software that automates simulations of ineslatic neutron scattering starting from ab-initio and ab-initio based calculations.


### DCS-Flow Requirements

* [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
* [DFTB+](https://www.dftbplus.org/)
* [Phonopy](https://phonopy.github.io/phonopy/)
* [Oclimax](https://neutrons.ornl.gov/sites/default/files/2018-NXS_Lecture_YQCheng_2.pdf)
   
### Linux Installation

* Install ASE:

```
pip install --upgrade --user ase
```

* Install DFTB+:

  1. [Download Slater-Koster files (parameters files for the DFTB method)](http://www.dftb.org/fileadmin/DFTB/public/slako-unpacked.tar.xz)

  2. [Download DFTB+ binaries](https://dftbplus.org/download/dftb-stable)

  3. Extract the file and put it in your home folder


* Installing Phonopy

```
pip install --upgrade --user phonopy
```

* Installing OCLIMAX

  1. OCLIMAX uses the DOCKER platform to run the application.
  [Please install it](https://www.docker.com/)

  2. Now download OCLIMAX

  ```  
  curl -sL https://sites.google.com/site/ornliceman/getoclimax | bash
  oclimax pull
  ```

* Installing DCS-Flow

```
git clone https://gitlab.com/lucassamir1/DCS-Flow.git
```

* Installing ChIMES

  If you want to use the machine learning method called ChIMES, please contact Nir Goldman (goldman14@llnl.gov) in order to receive the following files:

  1. DFTB+_ChIMES. This is a version of the DFTB+ package that has an interface to ChIMES. It should substitute the standart DFTB+ package.

  2. ChIMES binaries (chimes_lsq, chimes_md, lsq2.py). Put these files inside the git directory DCS-Flow/dcs/.

### Set environment variables

Add these lines to your configuration file (.bashrc). The following code uses example paths and must be edited according to your system.

```
export DFTB_PREFIX=/my_disk/my_name/slako/mio/mio-1-1/ # (path to Slako files)
export ASE_DFTB_COMMAND=/my_disk/my_name/dftbplus-20.1/bin/dftb+ > PREFIX.out # (path to dftb+)
export PATH=/my_disk/my_name/dftbplus-20.1/bin:$PATH (path to dftb+ files)
export PYTHONPATH=/my_disk/my_name/DCS-Flow:$PYTHONPATH #(path to DCS file)
export PATH=/my_disk/my_name/DCS-Flow/dcs:$PATH #(path to file within DCS folder)
```

### Usage

* Workflow

The workflow will relax the structure, create supercell displacements, calculate forces, run oclimax. The easiest way to start is through the command line interface

  1. Get the parameters file
   
  ```
  dcs workflow --get-params
  ```

  2. Define parameters in the file

  3. Run workflow

  ```
  dcs workflow
  ```


