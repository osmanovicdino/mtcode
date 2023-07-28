#!/bin/bash

# this code should be run from the base directory of the git

directory_path=$1
# param_file=$2 this file should be in the directory
g++ -fopenmp -std=c++17  ~/Chemistry/Code/mainNanotubeElasticShell_import.cpp -o ${directory_path}/a #copy to the directory we want 
fileparam=./params.csv
posfile=$(ls -t ${directory_path}/pos* | head -n 1)
orifile=$(ls -t ${directory_path}/ori* | head -n 1)
indfile=$(ls -t ${directory_path}/div* | head -n 1)
pfile=~/Chemistry/Code/IsocohedronI.csv
cp $pfile ${directory_path}
echo $fileparam
echo $posfile
echo $orifile
echo $pfile
cd $directory_path
./a $fileparam $posfile $orifile $indfile $pfile >log
