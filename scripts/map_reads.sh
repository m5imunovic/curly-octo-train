#!/bin/bash

# mamba install -c bioconda sra-tools minimap2 seqkit

# SRR_ENTRY="SRR9087597"
SRR_ENTRY="SRX5633451"

mkdir sra
# unpack sra file
mkdir fasta
mkdir fastq
# output sra entry to fasta directory using non-wrapped reads (i.e. 1-line == entire read)
fastq-dump -O fasta/ --fasta 0 sra/"$SRR_ENTRY"


mkdir paf
REF_DIR=references/chm13_v2.0
minimap2 -x map-hifi -H -t 20 -d "$REF_DIR"/ref.mmi "$REF_DIR"/chm13v2.0.fa.gz
# Fasta
# minimap2 -x map-hifi -H -c -t 20 "$REF_DIR"/ref.mmi fasta/"$SRR_ENTRY" > paf/alignment_"$SRR_ENTRY".paf
# Fastq
minimap2 -x map-hifi -H -c -t 20 "$REF_DIR"/ref.mmi fastq/"$SRR_ENTRY".fastq > paf/alignment_"$SRR_ENTRY".paf

