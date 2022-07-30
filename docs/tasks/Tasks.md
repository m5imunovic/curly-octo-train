This documents contains a list of tasks to implement specific steps of each pipeline
Each step should be specified with Heading 2 font. After finishing specific task it
should be linked to PR that implements it. The tasks implementing feature that is not
step specific should go under the `Generic Tasks` section. Use `Heading 3` for task
name

## Repo maintenance
### Improve handling of the OmegaConf objects
Modify the types in the calls of the functions that use configuration
to allow for the use of OmegaConf objects.


## Generic Tasks
~~### Create utility module for getting the commonly used paths in project
Use this to avoid hardcoded absolute paths in project~~

`#PR3:feature/create_util_path_handler`

### Create a module for the entire pipeline
It should contain all implemented steps in the pipeline.
Use top level config to control the execution.
Modify the simulation parameters to get the assembler graphs.

### Handle multiple chromosomes, instances of reads and assemblies
We should be able to handle multiple simulations of a chromosomes
in pipeline. Decide on the policy of reacting to existing data
(overwriting, skipping, aborting). 

### Improve logging of the pipeline
Use log instead of printing to the console.
Log the experiment config parameters (use hydra capabilities)

### Improve config integrity
Override the config values for the upcoming steps (e.g. species name)
if executed from the top level pipeline.


## Reference data 
~~### Create reference data simulator
Save simulated reference data into FASTA file(s)
Control via config file:
- chromosomes (name-length dict)
- GC content (float)
- name of the simulated species~~

`#PR4:feature/reference_simulator`

## Simulation
~~### Create simulation module
Based on PbSim2
Installation via script or python
Control via config file:
 - installation method
 - simulation parameters~~
 
 `#PR5:feature/read_simulator`

### Improve test
Add more tests
Mock complicated functions

## Assembly
~~### Create assembly module
Based on LJA
Installation via script or python
Control via config file:
 - installation method
 - assembly parameters~~

 `#PR5:feature/la_jolla_asm`

 ### Improve assembler testing
 -- Use fixtures to increase parameter testing coverage

## Ground truth data

## Loading Assembly Graph

## Graph feature extraction

