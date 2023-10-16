#!/bin/bash

# this code should be run from the base directory of the git

directory_path=$1
# param_file=$2 this file should be in the directory
g++ -fopenmp -std=c++17  ~/Chemistry/Code/mainNanotubeElasticShellPulsedImport.cpp -o ${directory_path}/a #copy to the directory we want 
fileparam=./param.csv
#cp $pfile ${directory_path}
echo $fileparam
# echo $posfile
# echo $orifile
# echo $indfile
# echo $pfile
cd $directory_path
export OMP_NUM_THREADS=8
par=`ls -l pos* | wc -l`
./a $fileparam $par >log
