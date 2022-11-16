#!/bin/bash
var=$1
gro=$2
cp $gro $var/sim.cpp
#g++ -g -rdynamic -std=c++17 -fopenmp $gro -o ${var}/angron
g++ -fopenmp -std=c++17 $gro -o ${var}/angron
