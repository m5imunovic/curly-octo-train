# Generating training datasets

In order to generate sequencing reads we are using two methods:

- simulation
- sampling

## Simulation

In simulation setting we are using already reconstructed reference to generate the data for ML experiment.
We utilize [PBSIM3](https://github.com/yukiteruono/pbsim3) to run the simulation. First we create the
simulation profile based on the reads that were used to assemble the reference. For each dataset subset
(train/val/test) we create different profiles. Additionally, for each sample, we use different seeds.
For each subset we internally also add the fixed offset to the seed:

- 10000 for `train` samples
- 20000 for `val` samples
- 30000 for `test` samples

## Sampling

In the sampling simulation we use real reads in order to generate the data for the ML experiment. First
we map the reads to chromosomes using [minimap](https://github.com/lh3/minimap2). For each chromosome we
perform subsampling, dropping certain percentage of reads, effectively reducing the coverage for each
particular chromosome. This way we can diversify data and model can learn in more challenging settings.

## Scenarios for building the datasets

Each dataset can be generated using the `scenario` file. Using scenario we can control how the samples are
created. For each sample we can control from which reference the chromosomes are bing used, the seed for
simulation, the (keep) probability for the sampling method. Each sample can be created from one or more
chromosomes.
