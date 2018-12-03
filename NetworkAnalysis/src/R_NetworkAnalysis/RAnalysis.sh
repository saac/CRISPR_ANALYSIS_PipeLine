#!/bin/bash


# Rscript SpacerAcquisitionNetwork.R Spacers-by-bacteria_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt
# Rscript ProtoSpacerAcquisitionNetwork.R Spacers-by-bacteria_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt
# Rscript ImmunityNetwork.R Bipartite_MATRIX_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt
# Rscript SpacerAcquisitionModularity.R Spacers-by-bacteria_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt

###################################


for i in $(ls ../Networks_$1/Spacers-by-bacteria_*)
do
        echo $i
#         Rscript ArgsTest.R $i 
        Rscript SpacerAcquisitionNetwork.R $i 
done

for i in $(ls ../Networks_$1/Protospacers-by-virus_*)
do
        echo $i
#         Rscript ArgsTest.R $i 
        Rscript ProtoSpacerAcquisitionNetwork.R $i

done

for i in $(ls ../Networks_$1/Bipartite_MATRIX_*)
do
        echo $i
#         Rscript ArgsTest.R $i 
        Rscript ImmunityNetwork.R $i

done

for i in $(ls ../Networks_$1/BipartieInfection_MATRIX_*)
do
        echo $i
#         Rscript ArgsTest.R $i 
        Rscript InfectionNetwork.R $i
done



# mkdir shuffled
for i in $(ls ../Networks_$1/Spacers-by-bacteria_*)
do
        echo $i
#         Rscript ArgsTest.R $i 
        Rscript SpacerAcquisitionModularity.R $i

done

# mkdir shuffled
for i in $(ls ../Networks_$1/Bipartite_MATRIX_*)
do
        echo $i
#         Rscript ArgsTest.R $i 
        Rscript ImmunityNetworkModularity.R $i

done


for i in $(ls ../Networks_$1/BipartieInfection_MATRIX_*)
do
        echo $i
#         Rscript ArgsTest.R $i 
        Rscript InfectionNetworkModularity.R $i
done


for i in $(ls ../Networks_$1/Protospacers-by-virus_*)
do
        echo $i
#         Rscript ArgsTest.R $i 
        Rscript ProtoSpacerAcquisitionModularity.R $i

done