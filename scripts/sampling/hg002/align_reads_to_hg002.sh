#!/bin/bash

READS_DIR="."

mkdir "$READS_DIR"/paf
for HG002_ENTRY in "HG002-rep1" "HG002-rep2" "HG002-rep3"
do
    minimap2 -x map-hifi -H -c -t 20 ./ref.mmi "$READS_DIR/fastq/$HG002_ENTRY.fastq.gz" > "$READS_DIR/paf/alignment_$HG002_ENTRY.paf"
done
