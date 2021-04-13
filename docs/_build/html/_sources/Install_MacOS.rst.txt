.. _installMac:

MacOS Installation
^^^^^^^^^^^^^^^^^^


* Install ASE:

.. code-block::

   pip3 install --upgrade --user ase


* 
  Install DFTB+:


  #. 
     `Download Slater-Koster files (parameters files for the DFTB method) <http://www.dftb.org/fileadmin/DFTB/public/slako-unpacked.tar.xz>`_

  #. 
     Download DFTB+

  .. code-block::

     git clone https://github.com/dftbplus/dftbplus.git
     cd dftbplus
     git submodule update --init --recursive

  _ Debugging Step: If getting an error related to WITH_OMP, open CMakeLists.txt file and add in the following line of code before if(WITH_OMP):

  .. code-block::

      option(WITH_OMP FALSE)


  #. Build DFTB+ (make sure to use your specific Fortran and C compilers)

  .. code-block::

     mkdir build
     cd build
     FC=gfortran CC=gcc cmake ..

  If configuration was successful

  .. code-block::

     ./dftb+
     make -j

  Test it

  .. code-block::

     ctest


  #. Install DFTB+

  .. code-block::

     make install


* Install Phonopy

.. code-block::

   pip3 install --upgrade --user phonopy


* 
  Install OCLIMAX


  #. 
     OCLIMAX uses the DOCKER platform to run the application.
     `Please install it <https://www.docker.com/>`_

  #. 
     Download OCLIMAX

  .. code-block::

     curl -sL https://sites.google.com/site/ornliceman/getoclimax | bash
     oclimax pull

* 
  Install DCS-Flow

.. code-block::

   git clone https://gitlab.com/lucassamir1/DCS-Flow.git


* 
  Installing ChIMES

  If you want to use the machine learning method called ChIMES, please contact Nir Goldman (goldman14@llnl.gov) in order to receive the following files:


  #. 
     DFTB+_ChIMES. This is a version of the DFTB+ package that has an interface to ChIMES. It should substitute the standart DFTB+ package.

  #. 
     ChIMES binaries (chimes_lsq, chimes_md, lsq2.py). Put these files inside the git directory DCS-Flow/dcs/.

Set environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^

Add these lines to your configuration file (~/.bash_profile). The following code uses example paths and must be edited according to your system.

.. code-block::

   export DFTB_PREFIX=/Users/my_name/slako/mio/mio-1-1/                                #(Path to Slako files)
   export ASE_DFTB_COMMAND=/Users/my_name/dftbplus/build/install/bin/dftb+ >PREFIX.out #(Path to dftb+)
   export PATH=/Users/my_name/dftbplus/build/install/bin:$PATH                         #(Path to dftb+)
   export PYTHONPATH=/Users/my_name/dftbplus/build/install/bin/dftb+:$PYTHONPATH       #(Python path to dftb+)
   export PATH=/Users/my_name/DCS-Flow/dcs:$PATH                               #(Path to DCS-Flow file)
   export PYTHONPATH=/Users/my_name/DCS-Flow:$PYTHONPATH                        #(Python path to DCS-Flow file)
   export PATH=/Users/my_name/.local/bin:$PATH                                         #(Path to ase file)
   export PYTHONPATH=/Users/my_name/.local/bin/ase:$PYTHONPATH                         #(Python path to ase file)

Usage
^^^^^


* Workflow

The workflow will relax the structure, create supercell displacements, calculate forces, run oclimax. The easiest way to start is through the command line interface


#. 
   Get the parameters file

   .. code-block::

      dcs workflow --get-params

#. 
   Define parameters in the file

#. 
   Run workflow

   .. code-block::

      dcs workflow
