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
    int runtime;
    if (argc == 6)
    {
        runtime = atof(argv[1]);
        packing_fraction = atof(argv[2]);
        int1 = atof(argv[3]);
        int2 = atof(argv[4]);
        int3 = atof(argv[5]);
    }
    else
    {
        runtime = 10000;
        packing_fraction = 0.01;
        int1 = 15.0;
        int2 = 10.0;
        int3 = 10.0;
    }
    matrix<double> params(28, 3);

    cout << packing_fraction << " " << int1 << " " << int2 << " " << int3 << endl;

    // for(int i = 0 ; i < 28 ; i++) {
    //     params(i, 0) =  10.0; //strength
    //     params(i, 1) =  1.4; //distance
    //     params(i, 2) = 0.927;
    // }

    int iter = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            params(iter, 0) = int1;
            params(iter, 1) = 1.4;
            params(iter, 2) = 0.927;
            iter++;
        }
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            if (j == 1)
            {
                params(iter, 0) = int2;
                params(iter, 1) = 1.4 * 0.75;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = int3;
                params(iter, 1) = 1.4 * 0.75;
                params(iter, 2) = 0.927;
            }
            iter++;
        }
    }

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {

            if (i == 0 && j == 0)
            {
                params(iter, 0) = int3;
                params(iter, 1) = 1.4 * 0.5;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1. * 0.5;
                params(iter, 2) = 0.927;
            }

            iter++;
        }
    }

    int n = 3000;
    int nt = 1000;
    TetrahedralWithBivalent c(params, nt, n);

    TetrahedralWithSingle c2(10.0, 1.4, 0.927, 10., 1.4, 0.927, 10.0, 1.4, 0.927, nt, n);

    TetrahedralPatch c3(20.0, 1.4, 0.927);

    //double packing_fraction = 0.02;

    double l = cbrt(pi * (double)n / (6. * packing_fraction));

    BindingModelBinary b(nt * 4);

    BindingModelSingle b2(0.998, 0.002);

    b.setup_equilibrium();

    Condensate A(l, n);

    // vector1<bool> pb(3, false);
    // cube geo(l, pb, 3);
    //A.obj->setgeometry(geo);

    A.setBindingModel(b2);

    A.setpots(c);

    A.setviscosity(0.1);

    double beta = 1.;

    A.obj->setkT(1. / beta);

    stringstream ss;
    ss << beta;

    string base = "_beta=";
    base += ss.str();

    double T;
    int TT;
    bool vv1, vv2, vv3;
    matrix<double> postemp = importcsv("/u/home/d/dinoo/Chemistry/Code/Basic/InitialConditions/evappos.csv", T, vv1);
    matrix<double> orienttemp = importcsv("/u/home/d/dinoo/Chemistry/Code/Basic/InitialConditions/evapori.csv", T, vv2);
    matrix<int> bindtemp = importcsv("/u/home/d/dinoo/Chemistry/Code/Basic/InitialConditions/evapbin.csv", TT, vv3);

    // matrix<double> postemp = importcsv("./Basic/InitialConditions/evappos.csv", T, vv1);
    // matrix<double> orienttemp = importcsv("./Basic/InitialConditions/evapori.csv", T, vv2);
    // matrix<int> bindtemp = importcsv("./Basic/InitialConditions/evapbin.csv", TT, vv3);

    matrix<double> newpostemp(n, 3);
    matrix<double> neworienttemp(n, 9);

    for(int i = 0  ; i < n ; i++) {
        newpostemp(i, 0) = postemp(i, 0);
        newpostemp(i, 1) = postemp(i, 1);
        newpostemp(i, 2) = postemp(i, 2);
        neworienttemp(i, 0) = orienttemp(i, 0);
        neworienttemp(i, 1) = orienttemp(i, 1);
        neworienttemp(i, 2) = orienttemp(i, 2);
        neworienttemp(i, 3) = orienttemp(i, 3);
        neworienttemp(i, 4) = orienttemp(i, 4);
        neworienttemp(i, 5) = orienttemp(i, 5);
        neworienttemp(i, 6) = orienttemp(i, 6);
        neworienttemp(i, 7) = orienttemp(i, 7);
        neworienttemp(i, 8) = orienttemp(i, 8);
    }

    A.obj->setdat(newpostemp);
    A.obj->setorientation(neworienttemp);
    BinaryBindStore bbs2;
    vector1<bool> iss(n);
    vector1<int> ist(n);
    for (int i = 0; i < n; i++)
    {
        iss[i] = (bool)bindtemp(0, i);
        ist[i] = bindtemp(1, i);
    }
    bbs2.isbound = iss;
    bbs2.boundto = ist;
    //Do processing to make sure everything is fine here

    A.run_singlebond_different_sizes_continue(runtime, 1000, nt, 0, bbs2, base);
    //A.run_singlebond(runtime, 1000, base);

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