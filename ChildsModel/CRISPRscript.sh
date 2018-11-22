#!/bin/bash

g++ -std=c++11 ChildsModel_8Clutser.cpp -o Crispr

max=$1

for i in $(seq 1 $max)
do
#     for j in 1 5 10 15 20
    for j in 1 10 20
    do
#         for k in 15 20
        for k in 15
        do
            sbatch ./CRISPR_B.sbatch $j 1e-7 5000 $k 10 $RANDOM
        done
    done
done

