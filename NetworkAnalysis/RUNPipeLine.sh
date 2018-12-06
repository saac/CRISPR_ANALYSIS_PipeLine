#!/bin/bash

# ./RUNPipeLine.sh T1 T2 FilePrefix-sed
# ./RUNPipeLine.sh 0 50 mu1e-7_initialDiffDp10_S10P15_R-1987

./FULLPlots.sh $1 $2


# mkdir Bacteria-abundance
# mkdir Phage-abundance

mkdir SpacerAcquisitionNetwork
mkdir ProtoSpacerAcquisitionNetwork
mkdir ImmunityNetwork

mkdir ImmunityNetwork-Modularity
mkdir SpacerAcquisitionModularity
mkdir ProtospacerAcquisitionModularity

mkdir InfectionNetwork
mkdir InfectionNetwork-Modularity


for i in $(seq $1 10 $2)
do
    echo $i
    ./PipeLine.sh $3 $i
    
    cd $i 
    
#     cp *Phage-abundance* ../Phage-abundance
#     cp *Bacteria-abundance* ../Bacteria-abundance
    
    cd R_NetworkAnalysis

    
    cp SpacerAcquisitionNetwork_$3_$i_*.png ../../SpacerAcquisitionNetwork
    cp ProtoSpacerAcquisitionNetwork_$3_$i_*.png ../../ProtoSpacerAcquisitionNetwork    
    cp ImmunityNetwork_$3_$i_*.png ../../ImmunityNetwork

    cp ImmunityNetwork-Modularity_$3_$i_*.png ../../ImmunityNetwork-Modularity    
    cp SpacerAcquisitionModularity_$3_$i_*.png ../../SpacerAcquisitionModularity
    cp ProtospacerAcquisitionModularity_$3_$i_*.png ../../ProtospacerAcquisitionModularity
    
    cp InfectionNetwork_$3_$i_*.png ../../InfectionNetwork 
    cp InfectionNetwork-Modularity_$3_$i_*.png ../../InfectionNetwork-Modularity      
    
    cd ../..
    rm -rf $i

done


mkdir $3_RESULTS_$1_$2

# mv Bacteria-abundance $3_RESULTS_$1_$2/
# mv Phage-abundance $3_RESULTS_$1_$2/
rm tmp
mv *_$3_Bacteria-abundance.png $3_RESULTS_$1_$2/
mv *_$3_Phage-abundance.png $3_RESULTS_$1_$2/

mv SpacerAcquisitionNetwork $3_RESULTS_$1_$2/
mv ProtoSpacerAcquisitionNetwork $3_RESULTS_$1_$2/
mv ImmunityNetwork $3_RESULTS_$1_$2/

mv ImmunityNetwork-Modularity $3_RESULTS_$1_$2/
mv SpacerAcquisitionModularity $3_RESULTS_$1_$2/
mv ProtospacerAcquisitionModularity $3_RESULTS_$1_$2/

mv InfectionNetwork $3_RESULTS_$1_$2/
mv InfectionNetwork-Modularity $3_RESULTS_$1_$2/



