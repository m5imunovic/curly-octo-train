#!/bin/bash

READS_ROOT="."
mkdir -p "$READS_ROOT"/index/

for SRR_ENTRY in "SRX7897685" "SRX7897686" "SRX7897687" "SRX7897688"
do
    for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 20 21 22 X
	do
	    echo "Processing $CHR..."
	    TARGET_ID=chr"$CHR"
	    OUTPUT_DIR=mapped/chr"$CHR"/$SRR_ENTRY/
	    mkdir -p "$OUTPUT_DIR"

	    INDEX_DIR="$READS_ROOT/index/$TARGET_ID"
	    mkdir -p "$INDEX_DIR"
	    # assumes tab delimited otherwise cut -d'<delimiter>'
	    grep "$TARGET_ID\s" "$READS_ROOT"/paf/alignment_"$SRR_ENTRY".paf | cut -f 1 > "$INDEX_DIR"/"$SRR_ENTRY".idx

	    # -f index file containing sequence ids - used as pattern
	    # -o output file name
	    # -w 0 - do not wrap

	    seqkit grep -f "$INDEX_DIR"/"$SRR_ENTRY".idx "$READS_ROOT"/fastq/"$SRR_ENTRY".fastq -o $OUTPUT_DIR/"$TARGET_ID".fastq -w 0
	done
done
