# CHM13 project

The [CHM13 project](https://github.com/marbl/CHM13) page contains detailed information about the
project and sequencing technologies that were used to assemble the reference.

## Simulation

In order to prepare the data for simulation we need a reference, reads and simulation profile. In
addition, make sure that the `vendor` tools are installed as describe in [install](TODO) instructions.

TODO: explain simulation procedure

## Sampling

If not already, make sure that the following prerequisite binaries are installed in environment:

```Bash
mamba install -c bioconda sra-tools minimap2 seqkit
```

In order to prepare the data for sampling we can use a set of scripts. In the instructions we are assuming
that the project data is located under `$data_dir` directory (default configuration is '/data').
Also, the assumption is that there is directory `references`, with already downloaded `chm13` reference
(`chm13v2.0.fa.gz` file) as explained in Simulation chapter.

First create a scratch folder and copy the scripts there from this repo:

```Bash
mkdir /data/scratch && cd /data/scratch
cp $PROJECT_ROOT/scripts/sampling/chm13/*.sh .
```

We need to download the reads first. This creates `fastq` directory with file blobs containing reads.

```Bash
bash download_chm13_reads.sh
```

Next thing we do is creating a reference mapping index, using recommended settings for hifi reads.
The output is index stored in `ref.mmi` file which we use in a subsequent step.

```Bash
bash index_chm13.sh
```

Once we have the index, we can align the reads to reference. This step produces alignments in `PAF`
format and saves them to local `paf`directory.

```Bash
bash align_reads_to_chm13.sh
```

Finally, we are able to produce the chromosome mapped reads. In the last step we are using alignments
to create index, so that we can map reads from blobs to respective chromosomes. These files are stored
in `mapped` directory, per chromosome.

```Bash
bash extract_map_chm13.sh
```

The contents should be copied to the `reads` directory:

```Bash
cp -r $data_dir/scratch/mapped $data_dir/references/homo_sapiens/chm13_v2/reads
```

The directory structure should look something like this:

```Bash
-> chm13_v2 tree reads/ | head -n 11
reads/
├── chr1
│   ├── SRX7897685
│   │   └── chr10.fastq
│   ├── SRX7897686
│   │   └── chr10.fastq
│   ├── SRX7897687
│   │   └── chr10.fastq
│   └── SRX7897688
│       └── chr10.fastq
├── chr2
...
```
