#!/bin/bash

var=$1
gro=$2
NN=`wc -l  $2 | awk '{ print $1 }'`
echo $NN

# for i in {1 .. $NN }
# do
#    echo $i
#    # echo `sed -n ${i}p ${gro} | expand -t 1`
# done

for (( c=1; c<=$NN+1; c++ ))
do
   echo `sed -n ${c}p ${gro} | expand -t 1`
   echo -n ""
   $var `sed -n ${c}p ${gro} | expand -t 1`
done



