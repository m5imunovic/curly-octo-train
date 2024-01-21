# curly-octo-train

## Install preprequisites

zlib (apt install zlib1g zlib1g-dev)
build-essentials (apt install build-essentials)

# Installation

We recommend the usage of mambaforge for setting up Python environment:

https://github.com/conda-forge/miniforge#mambaforge

```Bash
## Install conda
mamba init

## Create new environment (named menv)
mamba env create --file=environment.yaml
mamba activate menv
```

## Setup pre-commit hooks

```bash
pip install pre-commit
# Set up on the first usage
pre-commit install
pre-commit run --all-files
```

## Random reference genome step

In order to run the reference genome generator with default config execute the following command from command line:

```
PROJECT_ROOT="./" python src/reference.py
```

The reference config can be overridden with other configuration file:

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

## Genome assembly step

In order to run the assembly process on the reads generated by the genome sequencing step run the command:

```
PROJECT_ROOT="./" python src/assembly.py
```

The assembly is going to be created for all sequencing experiments of the species specified with `species_name` parameter. Currently it is expected that the reference and reads are per chromosome (i.e. one reference file for
chromosome and one reads file for that chromosome per sequencing experiment). The file is saved into the `assemblies/species_name/path_experiment_name` directory where path_experiment_name corresponds to the part of
the path of the reads file before the filename and after the root directory of the reads for that species.
One can limit the assembly to a specific set of paths using the cfg.asm.dir_filter config option as regex pattern.
All the files not matching the pattern will be filtered out. By default, all files are included, i.e. regex is left empty.

## Data annotation step

In order to create a machine learning set suitable for model experimentation it is necessary to export the
data files produced in assembly step into proper dataset format. This can be done by running the following command:

```
PROJECT_ROOT="./" python src/graph.py
```

In order to filter specific subset of assemblies that go into dataset we can use the `dir_filter` parameter:

```
PROJECT_ROOT="./" python src/graph.py  graph.set.dir_filter='S\[1-2\]'
```

## Evaluation of LJA

Branch `anton_development`, commit sha `9c81c179dcf8854981eaf843e2a44ba4c5607a7f` was used to perform evaluation of LJA metrics. In addition, following patch is applied
in order to add the export of initial GFA graph representation:

```diff
--- a/src/projects/error_correction/coverage_ec_stage.hpp
+++ b/src/projects/error_correction/coverage_ec_stage.hpp
@@ -32,6 +32,7 @@ CoverageEC(logging::Logger &logger, const std::experimental::filesystem::path &d
     io::SeqReader reader(reads_lib);
     readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
     printDot(dir / "initial_dbg.dot", Component(dbg));
+    printGFA(dir / "initial_dbg.gfa", Component(dbg), true);
     coverageStats(logger, dbg);
     if(debug) {
         PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, true);
diff --git a/src/projects/error_correction/topology_ec_stage.hpp b/src/projects/error_correction/topology_ec_stage.hpp
index 1709ec5..4c58f13 100644
--- a/src/projects/error_correction/topology_ec_stage.hpp
+++ b/src/projects/error_correction/topology_ec_stage.hpp
@@ -31,6 +31,7 @@ TopologyEC(logging::Logger &logger, const std::experimental::filesystem::path &d
     io::SeqReader reader(reads_lib);
     readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
     printDot(dir / "initial_dbg.dot", Component(dbg));
+    printGFA(dir / "initial_dbg.gfa", Component(dbg), true);
     if(debug) {
         DrawSplit(Component(dbg), dir / "before_figs", readStorage.labeler(), 25000);
         PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, false);
```
