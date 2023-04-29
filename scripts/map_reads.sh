#!/bin/bash

# mamba install -c bioconda sra-tools
# mamba install -c bioconda minimap2
# mamba install -c bioconda seq

SRR_ENTRY="SRR9087597"

mkdir fasta
mkdir sra
fastq-dump -O fasta/ --fasta 0 sra/"$SRR_ENTRY"
minimap2 -x map-hifi -H -t 20 -d ref.mmi references/chm13_v2.0/chm13v2.0.fa.gz
minimap2 -x map-hifi -H -c -t 20 ref.mmi fasta/"$SRR_ENTRY" > paf/alignment_"$SRR_ENTRY".paf

CHR=18
TARGET_ID=chr"$CHR"
# assumes tab delimited otherwise cut -d'<delimiter>'
mkdir paf
grep "$TARGET_ID" paf/"$SRR_ENTRY"_alignment.paf | cut -f 1 > paf/"$TARGET_ID"_reads.idx

# -f index file containing sequence ids - used as pattern
# -o output file name
# -w 0 - do not wrap

OUTPUT_DIR=reads/chr_"$CHR_real"/XX_XX/SXX/
mkdir -p $OUTPUT_DIR
seqkit grep -f paf/"$TARGET_ID"_reads.idx fasta/"$SRR_ENTRY".fasta -o $OUTPUT_DIR/"$TARGET_ID"_0001.fasta -w 0


