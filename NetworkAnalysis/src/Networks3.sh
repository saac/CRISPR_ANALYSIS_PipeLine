#!/bin/bash


# python BuildNetworks_v8.py mu1e-6_initialDiffDp10_S8P15_data-phage.txt mu1e-6_initialDiffDp10_S8P15_data-bact.txt T

# for i in {250,500,1000,1500,1750,2000,2500,3000}
# do
#         echo $i
#         python BuildNetworks_v11.py $1 $2 $i 
# done

python BuildNetworks_v12.py $1 $2 $4

mkdir Networks_$3
mv *Network*_R-* Networks_$3
mv *pacers-by* Networks_$3
mv Bipartite_MATRIX* Networks_$3
