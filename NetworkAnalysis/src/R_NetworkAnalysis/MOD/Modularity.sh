#!/bin/bash


echo $1

for i in $(ls $1)
do
#     echo $1/$i
    cp $1/$i .
done

for x in *.txt; do mv "$x" "${x%.txt}.sif"; done
#     mv *.txt *.sif

for i in $(ls *.sif)
do
    python Net2Infomap8.py $i tab nnI 1
done

mkdir MODULES

for i in $(ls *.net)
do
    ./Infomap --map --undirected -s$(echo $RANDOM) -N100 --tree --bftree --clu $i MODULES/
done

