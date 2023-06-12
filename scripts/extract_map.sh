#!/bin/bash

# SRR_ENTRY="SRR9087597"
SRR_ENTRY=SRX5633451

for CHR in 2 3 4 5 6 7 8 9 10 11 12 15 16 17 20 21 22 X
do
    # CHR=$1
    echo "Processing $CHR..."
    TARGET_ID=chr"$CHR"
    OUTPUT_DIR=reads/chr_"$CHR"_real/XX_XX/SXX/
    mkdir -p $OUTPUT_DIR

    # assumes tab delimited otherwise cut -d'<delimiter>'
    grep "$TARGET_ID\s" paf/alignment_"$SRR_ENTRY".paf | cut -f 1 > index/"$TARGET_ID".idx

    # -f index file containing sequence ids - used as pattern
    # -o output file name
    # -w 0 - do not wrap

    seqkit grep -f index/"$TARGET_ID".idx fastq/"$SRR_ENTRY".fastq -o $OUTPUT_DIR/"$TARGET_ID"_0001.fastq -w 0
done
