#!/bin/bash

READS_DIR="."

# SRR_ENTRIES:
#  SRX7897685
#  SRX7897686
#  SRX7897687
#  SRX7897688

# Fastq
mkdir "$READS_DIR"/paf
for SRR_ENTRY in "SRX7897685" "SRX7897686" "SRX7897687" "SRX7897688"
do
    minimap2 -x map-hifi -H -c -t 20 ./ref.mmi "$READS_DIR/fastq/$SRR_ENTRY.fastq" > "$READS_DIR/paf/alignment_$SRR_ENTRY.paf"
done
