#!/bin/bash 
# Assume that the fastq entry is in a single line
# paste prints four consecutive lines of fastq file in on row (tab delimited)
# cut takes only second column
# tr removes the new-lines
# wc counts characters
cat "$1" | paste - - - - | cut -f 2| tr -d '\n' | wc -c

# Seqk approach to count bases in fastq
seqkit sum -all "$1"
