#!/bin/bash
# get number of HPC compressed nucleotides in read (contained in read_1.txt file)
cat read_1.txt | tr -s 'A-Z' | tr -d '\n' | wc -m
