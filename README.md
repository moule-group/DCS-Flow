# Welcome to the repository for Adam Moule group - UC Davis

## Contents

* [Installation](https://gitlab.com/lucassamir1/adam-moule/-/tree/MacOSInstall#installation) 
* [Workflow](https://gitlab.com/lucassamir1/adam-moule/-/tree/MacOSInstall#workflow)
* [Documentation](https://gitlab.com/lucassamir1/adam-moule/-/tree/MacOSInstall#documentation)
* [Examples](https://gitlab.com/lucassamir1/adam-moule/-/tree/MacOSInstall#examples)


## Installation 

* [Linux Installation](https://gitlab.com/lucassamir1/adam-moule/-/blob/MacOSInstall/Install/Install_Linux.md)

* [MacOS Installation](https://gitlab.com/lucassamir1/adam-moule/-/blob/MacOSInstall/Install/Install_MacOS.md)


## Workflow 
- general workflow versus training workflow [e]

### Workflow Parameters Tips 
- add in guidance for how to choose krelax, fmax, geo, dim, kforce, mesh, task, e_unit values [e]


## Documentation
TODO:  summary of DCS purpose (caluclates phonon modes using correctd form of DFTB+ (general workflow), creates this corrected form by running DFT/md and machine learning process (training workflow)) Examples below show requirements for each part. [e]

DCS is a collection of the following scripts: 
* Relax: Optimizes structure
* Phonons: Calculates phonons
* Oclimax: Runs oclimax simulation
* Workflow: Workflow to relax structure, calculate phonons, and calculate the INS spectrum 

### Main Functions
#### Relax 
* ```relax(krelax=[6, 6, 6], fmax=0.05, geo=None, calc='dftbp'):```
    * Finds the geometry file and optimizes the structure using the specified calculator. 
    * The following input args are defined in the workflow parameters file:  
        * krelax (list, optional): Number of k points for relaxation. Defaults to [6, 6, 6].  
        * fmax (float, optional): Maximum allowed force for convergence between atoms. Defaults to 0.01.  
        * geo (str, optional): Geometry file or structure. 
            Allowed file types are .cif, .gen, .sdf, or .xyz. Defaults to None.  
        * calc (str, optional): Calculator used. Options are 'dftbp', 'chimes', or 'vasp'. Defaults to 'dftbp'.  

#### Train
* ```def train(dct=None):```  
    * Calls chosen calculator with a timer using specified parameters, else with default parameters.
    * The following input is required:
        * dct (dict, optional): Specified parameters for relax, md, and chimes functions. Defaults to None.

#### md
*  ```md(optgeo=None, calc='vasp', T=300, md_size=[1,1,1], steps=5000, time_step=1, dump_interval=100):```
    * Runs md using vasp or castep; if other calculator specified, raises error.
    * The following inputs are defined in the ____[q] file: 
        * optgeo (NoneType, optional): Optimized geometry file, only true if optgeo defined. Defaults to None.
        * calc (str, optional): Specifies calculator.Options are 'vasp' or 'castep'. Defaults to 'vasp'.
        * T (int, optional): Simulation temperature. Defaults to 300.
        * md_size (list, optional): Size of supercell. Defaults to [1,1,1].
        * steps (int, optional): Maximum number of ionic steps. Defaults to 5000.
        * time_step (int, optional): Md time step in fs. Defaults to 1. 
        * dump_interval (int, optional): Step size. Defaults to 100. 
    * Raises the following possible error [q] do I need to mention this? 
        * NotImplementedError: If calculator other than vasp or castep specified.

#### Phonons
* ```phonons(dim=[4, 4, 4], kforce=[1, 1, 1], mesh=[8, 8, 8], calc='dftbp'):```
    * Runs phonon supercell displacement calculations, populates 2-phonons file with results. 
    * The following input args are defined in the workflow parameters file: 
        * dim (list, optional): Dimensions of the supercell. Defaults to [4, 4, 4].
        * kforce (list, optional): Number of k points for force calculations. Defaults to [1, 1, 1].
        * mesh (list, optional): Uniform meshes for each axis. Defaults to [8, 8, 8].
        * calc (str, optional): Calculator used for task. Options are 'dftbp', 'chimes', 'vasp', or 'castep'. Defaults to 'dftbp'.
    * Raises the following possible error: 
        * NotImplementedError: Raised if 'calc' specified is not available. 

#### Oclimax
* ```oclimax(params=None, task=0, e_unit=0):```
    * Creates folder 3-oclimax within wd and write params file using default values if no dict exist in folder
    * The following input args are defined in the workflow parameters file: 
        * params (str, optional): Oclimax parameters file defined in write_params function. Defaults to None.
        * task (int, optional): Defines approximation method. 
            0:inc approx. 1:coh+inc. 2:single-xtal Q-E. 3:single-xtal Q-Q. Defaults to 0.
        * e_unit (int, optional): Defines energy unit. Defaults to 0.

#### Workflow
* ```workflow(dct=None):```
    * Calls all workflow functions with a timer using specified parameters, else with default parameters.
    * The following input dictionary is the workflow parameters files: 
        * dct (dict, optional): Specified parameters for relax, phonons, and oclimax functions. Defaults to None.

#### Chimes
* ```workflow(dct=None):```


## Examples

[q] so should I do on the super computer, personal computer too for DFTB?
[q] should i specifically say NERSC was used- is the run.py file special to NERSC or more general? 
[q] need to edit tesne- should it be instructing, or saying this is what was done 
---

__DFTB+ and Main Workflow for TCNQ on NERSC__:  

The following example shows the primary workflow using only dftb+ as the calculator. It was done using a super computer, but could also be run on a personal computer.  

First, a folder was created containing the geometry file (.cif, .gen, .sdf, or .xyz) and a run_tcnq.py file (for NERSC). This folder, named TCNQ, can be downloaded here: [q]  

The folder was uploaded to the super computer using the file transfer software Globus. [q] need to know how specific I should be - as opposed to "uploaded to NERSC, we used globus..."  

The super computer was accesed via the terminal and logged onto. Inside the TCNQ folder, the CNSS module was loaded, as shown.  

``` python
cd $SCRATCH 
cd TCNQ
module use /global/common/software/m2734/cnss/modulefiles
module load cnss
```  
Then the workflow parameters file, ```workflow_params.json```, was created using the following commands.  

``` python
cnss workflow --get-params
```  

The workflow parameters file was edited to match the following values. Check the Workflow Parameters Tips [e] [add link] for information on how to select values.  

```
{
    "krelax": [
        4,
        4,
        2
    ],
    "fmax": 0.05,
    "geo": null,
    "calc": "dftbp",
    "dim": [
        2,
        2,
        1
    ],
    "kforce": [
        1,
        1,
        1
    ],
    "mesh": [
        8,
        8,
        8
    ],
    "params": null,
    "task": 0,
    "e_unit": 0
}
```  

The TCNQ folder, or current directory, now has the structure file (tcnq.cif), the edited parameters file (workflow_params.json), and the run script (run_tcnq.py.). The run_tcnq.py file contains information for the NERSC super computer. The final lines contain the commands to be evaluated, in this case ```eval $'cnss -T workflow'```.  

The job is submitted and it's progress checked using the following commands:

``` python
sbatch run_tcnq.py
sqs
```  

Once the job has completed, the following files can be found in the TCNQ folder.  
```
1-optimization		3-oclimax		err.out			run_tcnq.py
2-phonons		PREFIX.out		out.out			workflow_params.json
```  

[q] not sure how to view the INS spectra with docker 

--- 

__DFTB+ and Main Workflow for TCNQ on PC__:  

The following example shows the primary workflow using only dftb+ as the calculator run on a personal terminal (as opposed to a super computer).  

First, a folder was created containing the geometry file (.cif, .gen, .sdf, or .xyz). This folder, named TCNQ, can be downloaded here: [q]   

In the TCNQ folder, the workflow parameters file, ```workflow_params.json```, was created using the following command.  

``` python
cnss workflow --get-params
```  

The workflow parameters file was edited to match the following values. Check the Workflow Parameters Tips [e] [add link] for information on how to select values.  

```
{
    "krelax": [
        4,
        4,
        2
    ],
    "fmax": 0.05,
    "geo": null,
    "calc": "dftbp",
    "dim": [
        2,
        2,
        1
    ],
    "kforce": [
        1,
        1,
        1
    ],
    "mesh": [
        8,
        8,
        8
    ],
    "params": null,
    "task": 0,
    "e_unit": 0
}
```  

The TCNQ folder, or current directory, now has the structure file (tcnq.cif) and the edited parameters file (workflow_params.json). The command to 

``` python
cnss -T workflow
```  

Once the job has completed, the following files can be found in the TCNQ folder.  
```
1-optimization		3-oclimax		err.out			run_tcnq.py
2-phonons		PREFIX.out		out.out			workflow_params.json
```  
[e] add how to view INS 

---

__Training and Main Workflow for TCNQ-TTF__:  

