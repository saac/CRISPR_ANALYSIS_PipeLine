#!/bin/bash


# for i in {150..250..10}
for i in {150..250..10}
do
    echo $i
#    ./SUPERMaster.sh mu1e-7_initialDiffDp1_S10P15_R-12499 $i
    ./SUPERMaster.sh mu1e-7_initialDiffDp1_S10P15_R-12499 $i
done


# for i in {100,200,500,1000,1900,2300,1750,2900,3500}
# do
#     echo $i
#     ./SUPERMaster.sh mu1e-7_initialDiffDp10_S10P15_R-1987 $i
# done