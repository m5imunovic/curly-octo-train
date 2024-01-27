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
