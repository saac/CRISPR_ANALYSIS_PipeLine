#!/bin/bash


# mkdir shuffled

for i in $(ls ../Networks_$1/Spacers-by-bacteria_*)
do
        Rscript SpacerAcquisitionAnalysis.R $i 
done


for i in $(ls ../Networks_$1/Protospacers-by-virus_*)
do
        Rscript ProtoSpacerAcquisitionAnalysis.R $i 
done


for i in $(ls ../Networks_$1/Bipartite_MATRIX_*)
do
        Rscript ImmunityNetworkAnalysis.R $i   
done


for i in $(ls ../Networks_$1/BipartieInfection_MATRIX_*)
do
        Rscript InfectionNetworkAnalysis.R $i        
done
