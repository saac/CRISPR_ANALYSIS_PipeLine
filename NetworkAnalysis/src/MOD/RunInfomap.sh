#!/bin/bash

# echo "${1%.txt}.sif"

cp $1 "${1%.txt}.sif"

python Net2Infomap8.py "${1%.txt}.sif" tab nnI 1

mkdir MODULES

./Infomap --map --undirected -s$(echo $RANDOM) -N100 --tree --bftree --clu "${1%.txt}.net" MODULES/

rm "${1%.txt}.pickle" 
# mv "${1%.txt}.pickle" 