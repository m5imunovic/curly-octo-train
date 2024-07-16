#!/bin/bash
READS_DIR="./fastq"
mkdir $READS_DIR

wget -O HG002-rep1.fastq.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifirevio/m84005_220827_014912_s1.hifi_reads.fastq.gz
mv HG002-rep1.fastq.gz $READS_DIR
wget -O HG002-rep2.fastq.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifirevio/m84011_220902_175841_s1.hifi_reads.fastq.gz
mv HG002-rep2.fastq.gz $READS_DIR
wget -O HG002-rep3.fastq.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifirevio/m84039_230418_213342_s3.hifi_reads.default.fastq.gz
mv HG002-rep3.fastq.gz $READS_DIR
