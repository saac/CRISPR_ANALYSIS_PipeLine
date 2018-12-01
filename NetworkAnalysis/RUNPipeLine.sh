#!/bin/bash

# ./RUNPipeLine.sh T1 T2 FilePrefix-sed
# ./RUNPipeLine.sh 0 50 mu1e-7_initialDiffDp10_S10P15_R-1987

./FULLPlots.sh 0 $2


# T1=$(($1 - 50))
# T2=$(($1 + 50))
# T1=$(($1 - 10))
# T2=$(($1 + 10))


# mkdir Bacteria-abundance
# mkdir Phage-abundance

mkdir SpacerAcquisitionNetwork
mkdir ProtoSpacerAcquisitionNetwork
mkdir ImmunityNetwork

mkdir ImmunityNetwork-Modularity
mkdir SpacerAcquisitionModularity
mkdir ProtospacerAcquisitionModularity


# for i in {1,2,3,4,5,6,7,8,9}
# do
#     echo $i
# #     ./SUPERMaster.sh mu1e-7_initialDiffDp1_S10P15_R-12499 $i
#     ./PipeLine.sh $3 $i
#     
# #     echo Networks_$2
#     cd $i 
#     
#     cp *Phage-abundance* ../Phage-abundance
#     cp *Bacteria-abundance* ../Bacteria-abundance
#     
#     cd R_NetworkAnalysis
# 
#     
#     cp SpacerAcquisitionNetwork_mu1e-7_$i_*.png ../../SpacerAcquisitionNetwork
#     cp ProtoSpacerAcquisitionNetwork_mu1e-7_$i_*.png ../../ProtoSpacerAcquisitionNetwork    
#     cp ImmunityNetwork_mu1e-7_$i_*.png ../../ImmunityNetwork
# 
#     cp ImmunityNetwork-Modularity_mu1e-7_$i_*.png ../../ImmunityNetwork-Modularity    
#     cp SpacerAcquisitionModularity_mu1e-7_$i_*.png ../../SpacerAcquisitionModularity
#     cp ProtospacerAcquisitionModularity_mu1e-7_$i_*.png ../../ProtospacerAcquisitionModularity
#     
#     cd ../..
# #     rm -rf $i
# 
# done



for i in $(seq $1 10 $2)
do
    echo $i
#     ./SUPERMaster.sh mu1e-7_initialDiffDp1_S10P15_R-12499 $i
    ./PipeLine.sh $3 $i
    
#     echo Networks_$2
    cd $i 
    
#     cp *Phage-abundance* ../Phage-abundance
#     cp *Bacteria-abundance* ../Bacteria-abundance
    
    cd R_NetworkAnalysis

    
    cp SpacerAcquisitionNetwork_mu1e-7_$i_*.png ../../SpacerAcquisitionNetwork
    cp ProtoSpacerAcquisitionNetwork_mu1e-7_$i_*.png ../../ProtoSpacerAcquisitionNetwork    
    cp ImmunityNetwork_mu1e-7_$i_*.png ../../ImmunityNetwork

    cp ImmunityNetwork-Modularity_mu1e-7_$i_*.png ../../ImmunityNetwork-Modularity    
    cp SpacerAcquisitionModularity_mu1e-7_$i_*.png ../../SpacerAcquisitionModularity
    cp ProtospacerAcquisitionModularity_mu1e-7_$i_*.png ../../ProtospacerAcquisitionModularity
    
    cd ../..
    rm -rf $i

done

# mkdir RESULTS_$1
# 
# mv Bacteria-abundance RESULTS_$1/
# mv Phage-abundance RESULTS_$1/
# 
# mv SpacerAcquisitionNetwork RESULTS_$1/
# mv ProtoSpacerAcquisitionNetwork RESULTS_$1/
# mv ImmunityNetwork RESULTS_$1/
# 
# mv ImmunityNetwork-Modularity RESULTS_$1/
# mv SpacerAcquisitionModularity RESULTS_$1/
# mv ProtospacerAcquisitionModularity RESULTS_$1/

# mkdir RESULTS_$1_$2
# 
# mv Bacteria-abundance RESULTS_$1_$2/
# mv Phage-abundance RESULTS_$1_$2/
# 
# mv SpacerAcquisitionNetwork RESULTS_$1_$2/
# mv ProtoSpacerAcquisitionNetwork RESULTS_$1_$2/
# mv ImmunityNetwork RESULTS_$1_$2/
# 
# mv ImmunityNetwork-Modularity RESULTS_$1_$2/
# mv SpacerAcquisitionModularity RESULTS_$1_$2/
# mv ProtospacerAcquisitionModularity RESULTS_$1_$2/


mkdir $3_RESULTS_$1_$2

mv Bacteria-abundance $3_RESULTS_$1_$2/
mv Phage-abundance $3_RESULTS_$1_$2/

mv SpacerAcquisitionNetwork $3_RESULTS_$1_$2/
mv ProtoSpacerAcquisitionNetwork $3_RESULTS_$1_$2/
mv ImmunityNetwork $3_RESULTS_$1_$2/

mv ImmunityNetwork-Modularity $3_RESULTS_$1_$2/
mv SpacerAcquisitionModularity $3_RESULTS_$1_$2/
mv ProtospacerAcquisitionModularity $3_RESULTS_$1_$2/

rm tmp
mv *_$3_Bacteria-abundance.png $3_RESULTS_$1_$2/
mv *_$3_Phage-abundance.png $3_RESULTS_$1_$2/



