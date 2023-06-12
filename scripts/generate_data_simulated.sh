#!/bin/bash
if [[ -z "$2" ]]; then
  if [[ $1 -eq "X" ]]; then
	SEED=$((23 + 100))
  else
	SEED=$(($1 + 100))
  fi
else
  SEED=$2
fi

echo "Using seed $SEED"
echo "Sequencing chromosome $1..."
PROJECT_ROOT="./" python src/sequencing.py --multirun seed=$SEED species_name@_global_=chr_$1_homunculus
echo "Assembling chromosome $1..."
PROJECT_ROOT="./" python src/assembly.py species_name@_global_=chr_$1_homunculus
echo "Creating dataset from chromosome $1..."
PROJECT_ROOT="./" python src/graph.py species_name@_global_=chr_$1_homunculus graph.set.name=test graph.debug=False
