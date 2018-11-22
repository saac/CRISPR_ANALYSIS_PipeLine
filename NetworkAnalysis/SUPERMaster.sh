#!/bin/bash


# ./SUPERMaster.sh "prefix_file" "Time"
# ./SUPERMaster.sh mu1e-7_initialDiffDp1_S10P15_R-12499 1711



mkdir $2
cp -r src/* $2
cp $1_*.txt $2

cd $2

# pwd

T1=$(($2 - 50))
T2=$(($2 + 50))

./FULLPlots.sh $T1 $T2
# pwd
./Master.sh $1 $2
# pwd

# ../$1_data-bact.txt $2
# ../$1_data-phage.txt $2

rm $1_*.txt

cd ..
# pwd
echo "FINISH!" $2

