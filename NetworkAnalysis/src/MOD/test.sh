#!/bin/bash


python BuildNetworks_v11.py mu1e-7_initialDiffDp1_S10P15_R-12499_data-phage.txt mu1e-7_initialDiffDp1_S10P15_R-12499_data-bact.txt $1
./RunInfomap.sh Similarity-Bacteria_Network_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_$1.txt 


