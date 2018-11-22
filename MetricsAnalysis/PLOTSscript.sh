#!/bin/bash


for tms in $(ls *time-series-data.txt)
do
	echo $tms 
        python TimeSeries1.py $tms
done


for i in $(ls *_data-phage.txt)
do 

    A="$(cut -d'_' -f1 <<<"$i")"
    B="$(cut -d'_' -f2 <<<"$i")"
    C="$(cut -d'_' -f3 <<<"$i")"
    D="$(cut -d'_' -f4 <<<"$i")"
    
    j=$(echo $A"_"$B"_"$C"_"$D)

    echo $i $j"_data-bact.txt"
    python DiversityIndices11.py $i $j"_data-bact.txt"
    
done

mkdir Plots
mv *.png Plots

mkdir RESULTS
mv *.txt RESULTS

./PLOTSscriptAVG2.sh
