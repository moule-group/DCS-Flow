# Welcome to the repository for Adam Moule group - UC Davis

## Contents

* Installation
* Workflow
* Documentation
* Examples 


## Installation 

* [Linux Installation](https://gitlab.com/lucassamir1/adam-moule/-/blob/MacOSInstall/Install/Install_Linux.md)

* [MacOS Installation](https://gitlab.com/lucassamir1/adam-moule/-/blob/MacOSInstall/Install/Install_MacOS.md)


## Workflow 


## Documentation
TODO:  summary of CNSS purpose (geom opt and create INS spectra)

CNSS is a collection of the following scripts: 
* Relax: Optimizes structure
* Phonons: Calculates phonons
* Oclimax: Runs oclimax simulation
* Workflow: Workflow to relax structure, calculate phonons, and calculate the INS spectrum 

### Functions 
(wouldn't necessarily include these, maybe for my thesis, but as an example )

 #### Relax 
 ```
 relax_done(fmax)
 ```
checks if relax.out file already exists (geom opt already done)
```
relax_structure(krelax, fmax, geo, mode)
```
uses specified mode (dftbp or vasp) to optmize geometry, returns if relax.out file already exists 
```
relax(krelax=[6, 6, 6], fmax=0.01, geo=None, calc='dftbp')
```
krelax: Number of k points for relaxation, e.g., 6 6 6 \n
fmax: Convergence criteria for forces \n
geo: Name of geometry file for structure (cif or gen extensions) \n 
calc: Calculator used; dftbp or vasp \n

THEN for example, i'd want to include the parameters, would they go elsewhere or under the func like i have 

        
#### Workflow

```write_params(task: int, e_unit: int)```
        writes out.params file
``` run_oclimax(params)```
        creates plot of Energy versus normalized intensity 
        params: JSON file with parameters for workflow





## Examples
