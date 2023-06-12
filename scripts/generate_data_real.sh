#!/bin/bash

for CHR in 1
#for CHR in 2 3 4 5 6 7 8 9
#for CHR in 10 11 12 13 14 15
#for CHR in 16 17 18 20 21 22 X
do
  echo "Assembling chromosome $CHR..."
  PROJECT_ROOT="./" python src/assembly.py species_name@_global_=chr_"$CHR"_real
  echo "Creating dataset from chromosome $CHR..."
  PROJECT_ROOT="./" python src/graph.py species_name@_global_=chr_"$CHR"_real graph.set.name=test graph.debug=False graph.set.dir_filter="chr$CHR/graph"
  echo "Eval La jolla"
  PROJECT_ROOT=./ python src/asm/eval_la_jolla.py asm_path=/trainingdata/exp/assemblies/chr_"$CHR"_real/04_30/S33/chr"$CHR"/ chr_n="$CHR"
done
