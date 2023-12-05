#!/bin/bash

# this code should be run from the base directory of the git

directory_path=$1
# param_file=$2 this file should be in the directory
#g++ -std=c++17  ~/Chemistry/Code/mainNanotubeElasticShell_import.cpp -o ${directory_path}/a #copy to the directory we want 
cp ~/Chemistry/Code/single_core_cont ${directory_path}/a
fileparam=./param.csv
posfile=$(ls ${directory_path}/pos* | tail -n 1)
orifile=./orient.csv
indfile=$(ls ${directory_path}/div* | tail -n 1)
pfile=./IsocohedronI.csv
#cp $pfile ${directory_path}
echo $fileparam
echo $posfile
echo $orifile
echo $indfile
echo $pfile
cd $directory_path
./a $fileparam $posfile $orifile $indfile $pfile >log
