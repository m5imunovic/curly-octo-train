#!/bin/bash

READS_ROOT="."
mkdir -p "$READS_ROOT"/index

for HG002_ENTRY in "HG002-rep1" "HG002-rep2" "HG002-rep3"
do
    for STRAND in "MATERNAL" "PATERNAL"
    do
        for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18 19 20 21 22 23
        do
            TARGET_ID=chr"$CHR"_"$STRAND"
            if [ $CHR == "23" ]; then
                if [ $STRAND == "MATERNAL" ]; then
                    TARGET_ID=chrX_"$STRAND"
                else
                    TARGET_ID=chrY_"$STRAND"
                fi
            fi
            echo "Processing $TARGET_ID..."
            OUTPUT_DIR=mapped/$TARGET_ID/$HG002_ENTRY/
            echo "mkdir -p $OUTPUT_DIR"
            mkdir -p "$OUTPUT_DIR"

            INDEX_DIR="$READS_ROOT/index/$TARGET_ID"
            mkdir -p "$INDEX_DIR"
            echo "mkdir -p $INDEX_DIR"
            # assumes tab delimited otherwise cut -d'<delimiter>'
            grep "$TARGET_ID\s" "$READS_ROOT"/paf/alignment_"$HG002_ENTRY".paf | cut -f 1 > "$INDEX_DIR"/"$HG002_ENTRY".idx

            # -f index file containing sequence ids - used as pattern
            # -o output file name
            # -w 0 - do not wrap

            echo "Extracting to $OUTPUT_DIR/$TARGET_ID.fastq"
            seqkit grep -f "$INDEX_DIR"/"$HG002_ENTRY".idx "$READS_ROOT"/fastq/"$HG002_ENTRY".fastq.gz -o $OUTPUT_DIR/"$TARGET_ID".fastq -w 0
        done
    done
done
