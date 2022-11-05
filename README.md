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


