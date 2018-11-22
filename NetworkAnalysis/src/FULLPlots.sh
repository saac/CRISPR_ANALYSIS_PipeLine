#!/bin/bash 

#./FULLPlots.sh T1 T2

# N=$1
D=1000
#D=$1

for i in $(ls *abundance.txt)
do
        N=$(awk  'BEGIN{max=0} NR>1 {if(($4)>max) max=($4)}END {print max}' $i)
        awk -v N="$N" -v D="$D" '{ if($4 < N/D) $3 = 0; print $0}' $i > tmp 
        python AREAplotFULL.py tmp $(tail -1 $i | awk '{ print $1 }') $i $1 $2
done
