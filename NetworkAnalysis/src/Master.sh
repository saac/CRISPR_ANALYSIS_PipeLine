#!/bin/bash

# ./Master.sh "prefix_file" "Time"
# ./Master.sh mu1e-7_initialDiffDp1_S10P15_R-12499 T


./Networks3.sh $1_data-phage.txt $1_data-bact.txt $1 $2

cd R_NetworkAnalysis/
./RAnalysis.sh $1
cd ..

# cd MOD/
# ./Modularity.sh ../Networks_$1/
# cd ..



