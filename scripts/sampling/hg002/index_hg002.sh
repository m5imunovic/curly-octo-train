#!/bin/bash

REF_DIR=$1
# -x preset for hifi reads
# -H homopolymer compress
# -t threads number
# -d dump to file
minimap2 -x map-hifi -H -t 20 -d ./ref.mmi "$REF_DIR"/hg002v1.0.1.fasta.gz
