# How to build apptainer

Optionally, export the APPTAINER_TMPDIR in order to use it as temporary working space when building.
Temporary space is also used when running containers in unprivileged mode, and performing some operations
on filesystems that do not fully support `--fakeroot`.

```Bash
export APPTAINER_TMPDIR=/scratch/apptainer/$USER
```

To build the `test.sandbox` from definition file use `dbgc.def`:

```Bash
apptainer build --sandbox --fix-perms --fakeroot test.sandbox dbgc.def
```

Optionally, one can modify the `test.sandbox` before freezing it into `sif` container file.

```Bash
apptainer shell --writable --fakeroot --nv --bind <local_path>:<remote_path>:<opts> test.sandbox/
```

Finally, create apptainer image:

```Bash
apptainer build curly.sif test.sandbox/
```

Do not forget to remove the `test.sandbox` in the end as it takes lots of storage space:

```Bash
rm -rf test.sandbox
```

The build assumes that there is a local copy of the repo in `$HOME/work/curly-octo-train`.
`PROJECT_ROOT` is located at `/work/curly-octo-train` in the container.
To run dataset creation execute:

```Bash
apptainer run --bind <ABSOLUTE_PATH_TO_HOME>:<ABSOULTE_PATH_TO_HOME> --nv curly.sif python3 src/experiment.py
```

The local `data` folder where the input files (references, simulation profiles) are typically mounted to this folder.

Currently, there is a dependency on scenarios being in project path and not in arbitrary config so
we need to mount that directory if we want to use non-commited scenario (this should be removed in
future versions once the scenario becomes standard hydra yaml configuration instead of JSON file).
One real example of apptainer command would then be:

```Bash
apptainer run \
        --bind $HOME/data:/data \
        --bind $HOME/data/cconfig/experiment/scenarios:/work/curly-octo-train/config/experiment/scenarios \
        $HOME/curly/curly.sif python3 src/experiment.py \
       --config-path /data/cconfig --config-name create_dataset_config.yaml \
        experiment.scenario.name=real_train_scenario
```
