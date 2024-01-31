# Pipeline specification

The pipeline for creating the machine learning dataset in Pytorch format consists of four steps:

1. Reference selection
2. Reads simulation
3. Assembly
4. Conversion to ML dataset
5. Maintenance and logging

The user specifies following parameters:

- Dataset name (new or existing dataset)
- Subset (train, val, test)
- Reference (simulated or supported real dataset)
- Reference subset (all chromosomes, subset of chromosomes, subset of each single chromosome)
- Profile(s) for simulation
- Simulation seed
- Assembler
- Features that should be used

For each sample we should keep the following info (dataset, reference selection, profile, features, creation date, git hash, config)

The pipeline should take care of scaling over the steps:

Reference step needs only one process, so it can take max number of available cores.
Read simulation also utilizes single thread, therefore, for each sample that we create schedule only a single process.

Assembly and graph stage are producing the most of the data. However, while the assembly step can utilize lot of cores,
the graph producing stage is mostly single threaded. Therefore it makes sense to split the computing power such that most
of the processors are dedicated to assembly. Once the first assembly is available we can start creating graph data sample
from it while a new assembly is being run in parallel.

The pipeline takes care of data cleanup and updating the dataset and metadata with logging. By default,
pipeline should remove all the files that are being created except when explicitly stated. As soon as the
graph stage for a sample is complete the cleanup should be scheduled as the pipeline creates a lot of data and
in case of large datasets this might cause it to break.

# Roles of each step

It is important to keep the responsibilities of each step isolated and independent as much as possible. This makes
configuration easier, and code simpler.
For example, the role of the "Reference" step is just to provide the existence of the genomic reference on the system.
It should not concern itself with the issues such as if this is going to be used once, multiple times, how the dataset
in which the data is flowing is produced or named, does the simulator only produces reads from subset or entire reference, etc.
