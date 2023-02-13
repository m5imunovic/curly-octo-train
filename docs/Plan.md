# Intro

This is a coarse plan of activities that needs to be achieved for reaching the first milestone.
The document should be updated whenever the certaing target objective is met or new objective
becomes visible.

## Pipeline implementation

Currently the following pipeline steps are required:

1. Setting up a reference

   - through random generator
   - from a known species

2. Creating the reads repository

   - through a simulator (e.g. PbSim2)
     - establish proper settings (coverage, error pct, etc)
     - create model for T2T dataset
   - using the sequencing data from a real sequencing effort

3. Running the JumboDBG step:

   - use reference sha or branch (development)
   - establish appropriate settings

4. Build the graph representation of the output data

   - Reconstruct graph (networkx or something similar)

5. Generate ground truth data

   - Based on multiplicity data

6. Run the graph analysis

   - Generate features (decide on ML library)

## Repository structure

1. Agree on PR process
2. Establish code quality checks
3. Maintain requirements.txt file
4. Set up CI
5. Maintain additional steps for the prerequisites of binaries in pipeline (e.g. PbSim2, JumboDBG)
