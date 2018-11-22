#!/bin/bash


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


for i in $(seq $1 10 $2)
do
    echo $i
#     ./SUPERMaster.sh mu1e-7_initialDiffDp1_S10P15_R-12499 $i
    ./SUPERMaster2.sh $3 $i
    
#     echo Networks_$2
    cd $i 
    
    cp *Phage-abundance* ../Phage-abundance
    cp *Bacteria-abundance* ../Bacteria-abundance
    
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


mkdir RESULTS_$1_$2

mv Bacteria-abundance RESULTS_$1_$2/
mv Phage-abundance RESULTS_$1_$2/

mv SpacerAcquisitionNetwork RESULTS_$1_$2/
mv ProtoSpacerAcquisitionNetwork RESULTS_$1_$2/
mv ImmunityNetwork RESULTS_$1_$2/

mv ImmunityNetwork-Modularity RESULTS_$1_$2/
mv SpacerAcquisitionModularity RESULTS_$1_$2/
mv ProtospacerAcquisitionModularity RESULTS_$1_$2/

