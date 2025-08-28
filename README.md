# mtcode
the CPU framework for the simulation of the microtubules.

the main file is mainMT.cpp, which provides a higher level control of the subsequent simulation. This file can be compiled using (for example):

g++ -std=c++11 mainMT.cpp -o a

where a is the name of the executable file created.

If one wishes to use CPU parallelization:

g++ -std=c++11 -fopenmp mainMT.cpp -o a

the output of the program is a csv of all the particle positions every 1000 timesteps of the simulation in whichever directory a is in.

The GPU method is found in the corresponding .cu files.

To run on GPU compile mainMT.cu using nvcc

01/27/2025

random_shuffle depreciated, replaced with shuffle and    

std::random_device rd; //to check whether this is too long
    std::mt19937 g(rd());

as the RNG. To be tested whether this decreases performance
