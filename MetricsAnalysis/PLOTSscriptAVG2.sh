#!/bin/bash


for j in 1 10 20
do
#   for k in 15 20
    for k in 15
    do

        awk 'FNR == 1 {next} {Bd[FNR]+=$5;Vd[FNR]+=$6;PDI[FNR]+=$9;T[FNR]++}END{for(i=2;i<=FNR;i++)print i-1,Bd[i]/T[i],Vd[i]/T[i],PDI[i]/T[i];}' RESULTS/$1"_initialDiffDp"$j"_S10"P"$k"_*time-series-data.txt > $1"_initialDiffDp"$j"_S10"P"$k"_AVERAGE_time-series-data.txt          
        python TimeSeriesAVG.py $1"_initialDiffDp"$j"_S10"P"$k"_AVERAGE_time-series-data.txt
        
        L=$(wc -l $1"_initialDiffDp"$j"_S10"P"$k"_AVERAGE_time-series-data.txt | awk '{print $1}')
        num=1
        L=$(($L+$num))
#         echo $L
        

        awk 'FNR == 1 {next} {Rch[FNR]+=$2;Sh[FNR]+=$3;if($4!="nan")Ev[FNR]+=$4;T[FNR]++}END{for(i=2;i<='$L';i++)print i-1,Rch[i]/T[i],Sh[i]/T[i],Ev[i]/T[i];}' RESULTS/$1"_initialDiffDp"$j"_S10"P"$k"*data-bact_DiversityIndices.txt > $1"_initialDiffDp"$j"_S10"P"$k"_AVERAGE_data-bact_DiversityIndices.txt
        awk 'FNR == 1 {next} {Rch[FNR]+=$2;Sh[FNR]+=$3;if($4!="nan")Ev[FNR]+=$4;T[FNR]++}END{for(i=2;i<='$L';i++)print i-1,Rch[i]/T[i],Sh[i]/T[i],Ev[i]/T[i];}' RESULTS/$1"_initialDiffDp"$j"_S10"P"$k"*data-phage_DiversityIndices.txt > $1"_initialDiffDp"$j"_S10"P"$k"_AVERAGE_data-phage_DiversityIndices.txt

        python DiversityIndicesAVG.py $1"_initialDiffDp"$j"_S10"P"$k"_AVERAGE_data-phage_DiversityIndices.txt $1"_initialDiffDp"$j"_S10"P"$k"_AVERAGE_data-bact_DiversityIndices.txt


    done
done

mkdir Averages
mv *AVERAGE* Averages