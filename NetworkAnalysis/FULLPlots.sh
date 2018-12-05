#!/bin/bash 

# N=$1
D=1000
#D=$1

for i in $(ls *abundance.txt)
do
        N=$(awk  'BEGIN{max=0} NR>1 {if(($4)>max) max=($4)}END {print max}' $i)  #max abundance of the file
        awk -v N="$N" -v D="$D" '{ if($4 < N/D) $3 = 0; print $0}' $i > tmp   #Temporary new file cut. Taking off the values less than N/D (D=1000)
        python AREAplotFULL.py tmp $(tail -1 $i | awk '{ print $1 }') $i $1 $2  #Run AREAplotFULL.py on the TeMPporary file and using as cutoff the last value of TIME in the original file
done
