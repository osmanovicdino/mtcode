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
#include <atomic>
#include <dirent.h>
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
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

int main(int argc, char **argv)
{

    srand(time(NULL));
    double packing_fraction;
    double int1;
    double int2;
    double int3;
    double beta;
    int runtime;
    if (argc == 4)
    {
        runtime = atof(argv[1]);
        packing_fraction = atof(argv[2]);
        beta = atof(argv[3]);
    }
    else
    {
        runtime = 10000;
        packing_fraction = 0.01;
        beta = 1.0;
    }

    cout << packing_fraction << " " << beta << endl;

    // for(int i = 0 ; i < 28 ; i++) {
    //     params(i, 0) =  10.0; //strength
    //     params(i, 1) =  1.4; //distance
    //     params(i, 2) = 0.927;
    // }

    int n = 9000;
    int nt = 1000;

    vector1<int> vec1(1);
    vec1[0] = 3;

    vector1<int> numb(1);

    numb[0] = nt;

    int tot = 0;
    for (int i = 0; i < 1; i++)
    {
        for (int j = i; j < 1; j++)
        {
            tot += vec1[i] * vec1[j];
        }
    }

    cout << tot << endl;

    matrix<double> params2(tot, 3);

    for (int i = 0; i < tot; i++)
    {
        params2(i, 0) = 10.0;
        params2(i, 1) = 1.4;
        params2(i, 2) = 1.1;
    }

    matrix<double> orient(3, 3);

    double nx1 = 0.0;
    double ny1 = 1.;
    double nz1 = 0.0;

    double nx2 = -sqrt(12. / 16.);
    double ny2 = -1. / 2.;
    double nz2 = 0.0;

    double nx3 = sqrt(12. / 16.);
    double ny3 = -1. / 2.;
    double nz3 = 0;

    matrix<double> asd(3, 3);

    asd(0, 0) = nx1;
    asd(0, 1) = ny1;
    asd(0, 2) = nz1;

    asd(1, 0) = nx2;
    asd(1, 1) = ny2;
    asd(1, 2) = nz2;

    asd(2, 0) = nx3;
    asd(2, 1) = ny3;
    asd(2, 2) = nz3;

    int iter2 = 0;
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < vec1[i]; j++)
        {
            orient(iter2, 0) = asd(j, 0);
            orient(iter2, 1) = asd(j, 1);
            orient(iter2, 2) = asd(j, 2);
            iter2++;
        }
    }

    GeneralPatch c4(vec1, numb, params2, orient);

    //double packing_fraction = 0.02;

    double l = cbrt(pi * (double)nt / (6. * packing_fraction));

    BindingModelBinary b(nt * 4);

    BindingModelSingle b2(0.998, 0.002);

    b.setup_equilibrium();

    Condensate A(l, nt);

    // vector1<bool> pb(3, false);
    // cube geo(l, pb, 3);
    //A.obj->setgeometry(geo);

    A.setBindingModel(b2);

    A.setpots(c4);

    A.setviscosity(0.1);

    //double beta = 1.;

    A.obj->setkT(1. / beta);

    stringstream ss;
    ss << beta;

    string base = "_beta=";
    base += ss.str();

    A.run_singlebond(runtime, 1000, base);

    //A.run_singlebond_different_sizes(runtime, 1000, nt, base);

    /*
int NN = 10000;
matrix<double> F(NN, 3);
matrix<double> T(NN, 3);

for(int i  = 0 ; i < 10000 ; i++) {
    F(i, 0) = (double)rand() / (double)(RAND_MAX);
    F(i, 1) = (double)rand() / (double)(RAND_MAX);
    F(i, 2) = (double)rand() / (double)(RAND_MAX);
    T(i, 0) = (double)rand() / (double)(RAND_MAX);
    T(i, 1) = (double)rand() / (double)(RAND_MAX);
    T(i, 2) = (double)rand() / (double)(RAND_MAX);
}

//int NN = A.obj->getN();

BinaryBindStore bbs;

int nh = (*A.pots).get_total_patches(NN);

vector1<bool> isbound(nh);

vector1<int> boundto(nh);

bbs.isbound = isbound;
bbs.boundto = boundto;
int ccc;
int num = 20;
matrix<int> boxes = A.obj->getgeo().generate_boxes_relationships(num, ccc);

matrix<int> *pairs = A.obj->calculatepairs(boxes, 3.5);

matrix<int> adj(28000,10);
vector1<int> len(28000,2);

for(int i = 0 ; i < 20000 ; i++) {
    //A.obj->CreateEdgeList(adj,len);
    //A.obj->advance_pos();
    cout << i << endl;
    A.obj->calculate_forces_and_torques3D_onlyone(*pairs, *(A.pots), bbs, *(A.bm), F, T);
    // stringstream ss;
    // ss<<i;
    // string res = "output"+ss.str();
    // outfunc(F,res);

 }
*/

    return 0;
}