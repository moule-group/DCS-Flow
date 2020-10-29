# Computational Neutron Scattering Simulation

CNSS is a software that automates simulations of ineslatic neutron scattering starting from ab-initio and ab-initio based calculations.


### CNSS Requirements

* [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
* [DFTB+](https://www.dftbplus.org/)
* [Phonopy](https://phonopy.github.io/phonopy/)
* [Oclimax](https://neutrons.ornl.gov/sites/default/files/2018-NXS_Lecture_YQCheng_2.pdf)

[Linux Installation](https://gitlab.com/lucassamir1/adam-moule/-/blob/MacOSInstall/cnss/README.md#linux-installation)

[MacOS Installation](https://gitlab.com/lucassamir1/adam-moule/-/blob/MacOSInstall/cnss/README.md#macos-installation)
   
### Linux Installation

* Install ASE:

```
pip install --upgrade --user ase
```

* Install DFTB+:

  1. [Download Slater-Koster files (parameters files for the DFTB method)](http://www.dftb.org/fileadmin/DFTB/public/slako-unpacked.tar.xz)

  2. Download DFTB+
  
  ```
  git clone https://github.com/dftbplus/dftbplus.git
  cd dftbplus
  git submodule update --init --recursive
  ```

  3. Build DFTB+ (make sure to use your specific Fortran and C compilers)

  ```
  mkdir build
  cd build
  FC=gfortran CC=gcc cmake ..
  ```

  If configuration was successful
  
  ```
  make -j
  ```

  Test it

  ```
  ctest -j2
  ```

  4. Install DFTB+

  ```
  make install
  ```

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

* Installing CNSS

```
git clone https://gitlab.com/lucassamir1/adam-moule.git
```

### Set environment variables

Add these lines to your configuration file (.bashrc for Linux, ~/.bash_profile for macOS). The following code uses example paths and must be edited according to your system.

```
export DFTB_PREFIX=/my_disk/my_name/slako/mio/mio-0-1/ # (path to Slako files)
export ASE_DFTB_COMMAND=/my_disk/my_name/dftbplus-20.1/bin/dftb+ > PREFIX.out # (path to dftb+)
export PATH=/my_disk/my_name/dftbplus-20.1/bin:$PATH (path to dftb+ files)
export PYTHONPATH=/my_disk/my_name/adam-moule/cnss:$PYTHONPATH #(path to CNSS file)
export PATH=/my_disk/my_name/adam-moule/cnss/cnss:$PATH #(path to file within CNSS folder)

```

### Usage

* Workflow

The workflow will relax the structure, create supercell displacements, calculate forces, run oclimax. The easiest way to start is through the command line interface

  1. Get the parameters file
   
  ```
  cnss workflow --get-params
  ```

  2. Define parameters in the file

  3. Run workflow

  ```
  cnss workflow
  ```


