# curly-octo-train

## Random reference genome step
In order to run the reference genome generator with default config execute the following command from command line:

```
PROJECT_ROOT="./" python src/reference.py 
```
The reference config can be overriden with other configuration file:
```
PROJECT_ROOT="./" python src/reference.py species_name@_global_=random_species_small_c1
```
If you want to create a new species config make sure to also provide the matching config in the `species_name` directory

You can also use `hydra --multirun` option to execute the multiple jobs with different configs simultaneously:
```
PROJECT_ROOT="./" python src/reference.py --multirun species_name@_global_=random_species_small_c5,random_species_medium_c1,random_species_medium_c5
```
In case that you are using some other sweeper where parallel jobs are supported it is recommended to adjust the `reads.threads` value accordingly.

## Simulate genome sequencing step

In order to simulate sequencing process with default config execute the following command from command line:

```
PROJECT_ROOT="./" python src/sequencing.py 
```
One can use the seed parameter in order to simulate multiple sequencing runs from the same reference genome:
```
PROJECT_ROOT="./" python src/sequencing.py --multirun seed=1,2,3,4,5 
```

