#!/bin/bash


for i in {200,500,750,1000,1250,1500,1750,2000,2500}
do
#     echo $i
    ./SUPERDuperMaster2.sh $i mu1e-7_initialDiffDp1_S10P15_R-12499
done
