Documentation
-------------

DCS-Flow is a collection of the following scripts: 


* Relax: Optimizes structure.
* Phonons: Calculates phonons modes with the supercell method.
* Oclimax: Runs oclimax simulation creating a INS sprectrum. 
* MD: Runs molecular dynamics simulations.
* Chimes: Creates the coefficients for the Chebyshev Interaction Model for Efficient Simulation. 
* Train: Automated workflow that calls the functions necessary to create the ChIMES coefficients.
* Workflow: Automates the main workflow functions to relax structure, calculate phonons, and calculate the INS spectrum. 

DCS-Flow has a command line interface implemented. Examples for how to use it are included under each main function. 

Main Functions
^^^^^^^^^^^^^^

Relax
~~~~~


* ``relax(krelax=[6, 6, 6], fmax=0.05, geo=None, calc='dftbp', T=5):``

  * Finds the geometry file and optimizes the structure using the specified calculator (Populates 1-optimization folder with results). 
  * The following input args are defined in the workflow parameters file:  

    * krelax (list, optional): Number of k points for relaxation. Defaults to [6, 6, 6].  
    * fmax (float, optional): Maximum allowed force for convergence between atoms. Defaults to 0.01.  
    * geo (str, optional): Geometry file or structure. 
        Allowed file types are .cif, .gen, .sdf, or .xyz. Defaults to None.  
    * calc (str, optional): Calculator used. Options are 'dftbp', 'chimes', 'castep', or 'vasp'. Defaults to 'dftbp'.
    * T (int, optional): Simulation temperature. Defaults to 5K. Only used for DFTB+ and ChIMES.

.. code-block::

   dcs relax --krelax 6 6 6 --fmax 0.05 --geo TCNQ.cif --calc dftbp --temp 5

Phonons
~~~~~~~


* ``phonons(dim=[4, 4, 4], kforce=[1, 1, 1], mesh=[8, 8, 8], calc='dftbp', T=5):``

  * Runs phonon supercell displacement calculations, populates 2-phonons folder with results. 
  * The following input args are defined in the workflow parameters file: 

    * dim (list, optional): Dimensions of the supercell. Defaults to [4, 4, 4].
    * kforce (list, optional): Number of k points for force calculations. Defaults to [1, 1, 1].
    * mesh (list, optional): Uniform meshes for each axis. Defaults to [8, 8, 8].
    * calc (str, optional): Calculator used for task. Options are 'dftbp', 'chimes', 'vasp', or 'castep'. Defaults to 'dftbp'.
    * T (int, optional): Simulation temperature. Defaults to 5K. Only used for DFTB+ and ChIMES.
.. code-block::

   dcs phonons --dim 4 4 4 --kforce 1 1 1 --mesh 8 8 8 --calc dftbp --temp 5

Oclimax
~~~~~~~


* ``oclimax(params=None, task=0, e_unit=0):``

  * Creates folder 3-oclimax, writes oclimax parameters file in folder and runs oclimax simulation.
  * The following input args are defined in the workflow parameters file: 

    * params (str, optional): Oclimax parameters file name if exists. Otherwise, it will be created in write_params function. Defaults to None.
    * task (int, optional): Defines approximation method. 
        0:inc approx. 1:coh+inc. 2:single-xtal Q-E. 3:single-xtal Q-Q. Defaults to 0.
    * e_unit (int, optional): Defines energy unit. Defaults to 0.

.. code-block::

   dcs oclimax --task 0 --e_unit 0

MD
~~


* ``md(optgeo=None, calc='vasp', T=300, md_size=[1,1,1], steps=5000, time_step=1, dump_interval=100):``

  * Runs molecular dynamics simulation using vasp or castep (Creates 2-molecular_dynamics folder inside 0-train).
  * The following inputs are defined in the training parameters file (train_params.json): 

    * optgeo (NoneType, optional): Optimized geometry file, only true if optgeo defined. Defaults to None.
    * calc (str, optional): Specifies calculator.Options are 'vasp' or 'castep'. Defaults to 'vasp'.
    * T (int, optional): Simulation temperature. Defaults to 300.
    * md_size (list, optional): Size of supercell. Defaults to [1,1,1].
    * steps (int, optional): Maximum number of ionic steps. Defaults to 5000.
    * time_step (int, optional): Md time step in fs. Defaults to 1. 
    * dump_interval (int, optional): Step size of frames to be saved in the trajectory file. Defaults to 100. 

.. code-block::

   dcs md --calc vasp --T 300 --md_size 1 1 1 --steps 5000 --time_step 1 --dump_interval 100

Chimes
~~~~~~


* ``chimes(trajfile=None, b2=12, b3=8, T=5):``

  * Calculates force difference between DFT and DFTB (training set), and fits the Chebyshev polynomials coefficients. Creates 3-chimes folder (inside 0-train directory) and writes params.txt file.
  * The following inputs are required:

    * trajfile (list, optional): Trajectory file output from md simulation. Defaults to None.
    * b2 (int, optional): Second body order of Chebyshev polynomial. Defaults to 12.
    * b3 (int, optional): Third body order of Chebyshev polynomial. Defaults to 8.
    * T (int, optional): Temperature for simulation in Kevin. Defaults to 5.

.. code-block::

   dcs chimes --b2 12 --b3 8 --T 5

Train
~~~~~


* ``train(dct=None):``  

  * Calls functions related to the training workflow (relax, md, chimes) with a timer using specified parameters in train_params.json, else with default parameters. Creates 0-train directory.
  * The following input is required:

    * dct (dict, optional): JSON file with specified parameters for relax, md, and chimes functions. Defaults to 'train_params.json'.

.. code-block::

   dcs train

Workflow
~~~~~~~~


* ``workflow(dct=None):``

  * Calls all workflow functions (relax, phonons, oclimax) with a timer using specified parameters in workflow_params.json, else with default parameters.
  * The following input dictionary is the workflow parameters files: 

    * dct (dict, optional): JSON file with specified parameters for relax, phonons, and oclimax functions. Defaults to 'workflow_params.json'.

.. code-block::

   dcs workflow


