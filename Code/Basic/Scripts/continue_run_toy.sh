#!/bin/bash

# this code should be run from the base directory of the git

directory_path=$1
# param_file=$2 this file should be in the directory
g++ -fopenmp -std=c++17  ~/Chemistry/Code/mainNanotubeElasticShellSF_import.cpp -o ${directory_path}/a #copy to the directory we want 
fileparam=./param.csv
posfile=$(ls -t ${directory_path}/pos* | head -n 1)
orifile=$(ls -t ${directory_path}/pols* | head -n 1)
indfile=$(ls -t ${directory_path}/div* | head -n 1)
#cp $pfile ${directory_path}
echo $fileparam
echo $posfile
echo $orifile
echo $indfile
echo $pfile
cd $directory_path
export OMP_NUM_THREADS=6
./a $fileparam $posfile $orifile $indfile >log
