This documents contains a list of tasks to implement specific steps of each pipeline
Each step should be specified with Heading 2 font. After finishing specific task it
should be linked to PR that implements it. The tasks implementing feature that is not
step specific should go under the `Generic Tasks` section. Use `Heading 3` for task
name

## Repo maintenance

## Generic Tasks
~~### Create utility module for getting the commonly used paths in project
Use this to avoid hardcoded absolute paths in project~~

`#PR3:feature/create_util_path_handler`

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

### Improve handling of the OmegaConf objects
Modify the types in the calls of the functions

## Assembly
### Create assembly module
Based on LJA
Installation via script or python
Control via config file:
 - installation method
 - assembly parameters

## Ground truth data

## Loading Assembly Graph

## Graph feature extraction

