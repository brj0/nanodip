#!/bin/bash

out=all_clf.csv

# Header
nm_groups=$( echo ,meth_grp{1..10} | sed 's/ //g' )
nm_probs=$( echo ,prob{1..10} | sed 's/ //g' )
echo case,clf,acc,time$nm_groups$nm_probs >> $out

clf=(RandomForestClassifier KNeighborsClassifier MLPClassifier)

for i in 0 1 2;
do
    for file in $@;
    do
        start=$((16*i))
        case_no=${file::-45}
        methyl_grps=$( cat $file | head -$(( start + 15 )) | tail -10 | awk '{print $1}' | paste -s -d ,)
        prob=$( cat $file | head -$(( start + 15 )) | tail -10 | awk '{print $3}' | paste -s -d ,)
        acc=$( cat $file | head -$(( start + 4 )) | tail -1 | awk '{print $3}' | paste -s -d ,)
        time=$( cat $file | head -$(( start + 2 )) | tail -1 | awk '{print $5}' | paste -s -d ,)
        echo $case_no,${clf[i]},$acc,$time,$methyl_grps,$prob >> $out
    done
done
