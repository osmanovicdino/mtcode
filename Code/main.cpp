#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/ioctl.h> 
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <random>
#include <mutex>
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "Basic/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"
//#include "MDBase/potential.h"
#include "MDBase/MD.h"
#include "MDBase/Langevin.h"
//#include "MDBase/LangevinR.h"
#include "Condensate/Condensate.h"

// #include "NCGasR.h"
// #include "Microtubule.h"



//#include "MDGPU.cu"

using namespace std;



int main(int argc, char** argv) {

srand (time(NULL));

if(argc!=5) error("incorrect number of arguments passed to main file");

double packing_fraction = atof(argv[1]);
double int1 = atof(argv[2]);
double int2 = atof(argv[3]);
double int3 = atof(argv[4]);

matrix<double> params(28,3);

cout << packing_fraction << " " << int1 << " " << int2 << " " << int3 << endl;

// for(int i = 0 ; i < 28 ; i++) {
//     params(i, 0) =  10.0; //strength 
//     params(i, 1) =  1.4; //distance
//     params(i, 2) = 0.927;
// }

int iter = 0;
for(int i = 0  ; i < 4 ; i++) {
    for(int j = 0  ; j < 4 ; j++) {
        params(iter,0) = int1;
        params(iter,1) = 1.4;
        params(iter,2) =  0.927;
        iter++;
    }
}

for(int i = 0  ; i < 4 ; i++) {
    for(int j = 0 ; j < 2 ; j++) {
            params(iter, 0) = int2;
            params(iter, 1) = 1.4;
            params(iter, 2) = 0.927;
        
        iter++;
    }
}

for (int i = 0; i < 2; i++)
{
    for (int j = 0; j < 2; j++)
    {

            params(iter, 0) = int3;
            params(iter, 1) = 1.4;
            params(iter, 2) = 0.927;
        

        iter++;
    }
}


int n = 10000;
int nt = 2000;
TetrahedralWithBivalent c(params,nt,n);

TetrahedralWithSingle c2(10.0, 1.4, 0.927, 10., 1.4, 0.927, 10.0, 1.4, 0.927, nt, n);


TetrahedralPatch c3(10.0,1.4,0.927);

//double packing_fraction = 0.02;

double l = cbrt(pi * (double)n / (6. * packing_fraction));

BindingModelBinary b(nt*4);

b.setup_equilibrium();

Condensate A(l, n);

A.setBindingModel(b);

A.setpots(c);

A.setviscosity(1.0);

double beta = 1.0;

A.obj->setkT(1. / beta);

stringstream ss;
ss << beta;

string base = "_beta=";
base += ss.str();

A.run_singlebond(10000, 1000, base);

return 0;
}