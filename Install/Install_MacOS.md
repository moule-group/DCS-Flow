# Davis Computational Spectroscopy workflow (DCS-Flow)

DCS-Flow is a software that automates simulations of ineslatic neutron scattering starting from ab-initio and ab-initio based calculations.


### DCS-Flow Requirements

* [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
* [DFTB+](https://www.dftbplus.org/)
* [Phonopy](https://phonopy.github.io/phonopy/)
* [Oclimax](https://neutrons.ornl.gov/sites/default/files/2018-NXS_Lecture_YQCheng_2.pdf)

### MacOS Installation

* Install ASE:

```
pip3 install --upgrade --user ase
```

* Install DFTB+:

  1. [Download Slater-Koster files (parameters files for the DFTB method)](http://www.dftb.org/fileadmin/DFTB/public/slako-unpacked.tar.xz)

  2. Download DFTB+
  
  ```
  git clone https://github.com/dftbplus/dftbplus.git
  cd dftbplus
  git submodule update --init --recursive
  ```
  _ Debugging Step: If getting an error related to WITH_OMP, open CMakeLists.txt file and add in the following line of code before if(WITH_OMP):
  ```
   option(WITH_OMP FALSE)
   ```

  3. Build DFTB+ (make sure to use your specific Fortran and C compilers)

  ```
  mkdir build
  cd build
  FC=gfortran CC=gcc cmake ..
  ```

  If configuration was successful
  
  ```
  ./dftb+
  make -j
  ```

  Test it

  ```
  ctest
  ```

  4. Install DFTB+

  ```
  make install
  ```
  

* Install Phonopy

```
pip3 install --upgrade --user phonopy
```

* Install OCLIMAX

  1. OCLIMAX uses the DOCKER platform to run the application.
  [Please install it](https://www.docker.com/)

  2. Download OCLIMAX

  ```  
  curl -sL https://sites.google.com/site/ornliceman/getoclimax | bash
  oclimax pull
  ```

* Install DCS-Flow

```
git clone https://gitlab.com/lucassamir1/DCS-Flow.git
```


### Set environment variables

Add these lines to your configuration file (.bashrc for Linux, ~/.bash_profile for macOS). The following code uses example paths and must be edited according to your system.

```
export DFTB_PREFIX=/Users/my_name/slako/mio/mio-1-1/                                #(Path to Slako files)
export ASE_DFTB_COMMAND=/Users/my_name/dftbplus/build/install/bin/dftb+ >PREFIX.out #(Path to dftb+)
export PATH=/Users/my_name/dftbplus/build/install/bin:$PATH                         #(Path to dftb+)
export PYTHONPATH=/Users/my_name/dftbplus/build/install/bin/dftb+:$PYTHONPATH       #(Python path to dftb+)
export PATH=/Users/my_name/DCS-Flow/dcs:$PATH                               #(Path to DCS-Flow file)
export PYTHONPATH=/Users/my_name/DCS-Flow:$PYTHONPATH                        #(Python path to DCS-Flow file)
export PATH=/Users/my_name/.local/bin:$PATH                                         #(Path to ase file)
export PYTHONPATH=/Users/my_name/.local/bin/ase:$PYTHONPATH                         #(Python path to ase file)

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


