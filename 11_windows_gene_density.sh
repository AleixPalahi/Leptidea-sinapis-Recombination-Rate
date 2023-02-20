#!/bin/bash

awk '
function ceil(x, y){
 y=int(x);
 return(x>y?y+1:y)};
 {window_size=1000000;
 a[$1 "\t" ceil($2/window_size)] += ($3-$2)}
END {
 for (wind in a) print wind "\t" a[wind]}' cds_sorted.bed | sort -k1,1 -k2,2n > cds_1Mb.txt





