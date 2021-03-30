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
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
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
    //signal(SIGSEGV, handler);

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
        runtime = 100000;
        packing_fraction = 0.05;
        int1 = 12.0;
        int2 = 22.0;
        int3 = 7.0;
    }

    cout << packing_fraction << " " << int1 << " " << int2 << "  " << int3 << endl;

    // for(int i = 0 ; i < 28 ; i++) {
    //     params(i, 0) =  10.0; //strength
    //     params(i, 1) =  1.4; //distance
    //     params(i, 2) = 0.927;
    // }
    int m1 = 2000;
    int m2 = 6000;
    int n = 10000;

    // int m1  = 2;
    // int m2 = 6;
    // int n = 10;

    BindingModelTernary b(m1 * 4, m1*4 + 4 * (m2-m1));

    // b.setup(0.99,0.01,0.01,0.01,0.0,0.0,
    // 0.,
    // 0.0,
    // 0.0,
    // 0.0,
    // 0.0,
    // 0.0,
    // 0.0,
    // 0.0,
    // 0.0);

    double unbinding_rebar =  -2.0;

    b.setup(0.99, 0.01, 0.99, 0.2, 0.99, 0.2,
            0.,
            0.,
            0.,
            0.0,
            unbinding_rebar,
            0.0,
            0.0,
            0.0,
            0.0);

    // vector1<int> cc(4);

    // for(int k = 0; k < 100000 ; k++) {
    // bool after1, after2,after3;
    // b.triplet(false,false,true,true,true,true, 6000,20000,28000,after1,after2,after3);

    // if(after1) {
    //     cc[0]++;
    // }
    // else if(after2) {
    //     cc[1]++;
    // }
    // else if(after3) {
    //     cc[2]++;
    // }
    // else {
    //     cc[3]++;
    // }
    // }
    // cout << cc << endl;
    // pausel();

    
    // pausel();

    vector1<int> vec1(3);
    vec1[0] = 4;
    vec1[1] = 4;
    vec1[2] = 2;

    vector1<int> numb(3);

    numb[0] = m1;
    numb[1] = m2;
    numb[2] = n;

    int tot = 4*4+4*2+4*4+3*2+4*2+4*4;
    matrix<double> params(tot, 3);

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
        for (int j = 0; j < 4; j++)
        {

            if (i == 0 && j == 0)
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1.4;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1.4;
                params(iter, 2) = 0.927;
            }

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
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {

            if (i == 0 && j == 0)
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1 * 0.5;
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

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 2; j++)
        {

            if (i == 0 && j == 0)
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1.4 * 0.75;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = 0.0;
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
            if (j == 1)
            {
                params(iter, 0) = int3;
                params(iter, 1) = 1.4 * 0.5;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1.4 * 0.5;
                params(iter, 2) = 0.927;
            }
            iter++;
        }
    }

    matrix<double> orient(10,3);

    double nx1 = sqrt(8. / 9.);
    double ny1 = 0.;
    double nz1 = -1. / 3.;

    double nx2 = -sqrt(2. / 9.);
    double ny2 = sqrt(2. / 3.);
    double nz2 = -1. / 3.;

    double nx3 = -sqrt(2. / 9.);
    double ny3 = -sqrt(2. / 3.);
    double nz3 = -1. / 3.;

    double nx4 = 0;
    double ny4 = 0;
    double nz4 = 1.;



    orient(0, 0) = nx1;
    orient(0, 1) = ny1;
    orient(0, 2) = nz1;

    orient(1, 0) = nx2;
    orient(1, 1) = ny2;
    orient(1, 2) = nz2;

    orient(2, 0) = nx3;
    orient(2, 1) = ny3;
    orient(2, 2) = nz3;

    orient(3, 0) = nx4;
    orient(3, 1) = ny4;
    orient(3, 2) = nz4;

    orient(4, 0) = nx1;
    orient(4, 1) = ny1;
    orient(4, 2) = nz1;

    orient(5, 0) = nx2;
    orient(5, 1) = ny2;
    orient(5, 2) = nz2;

    orient(6, 0) = nx3;
    orient(6, 1) = ny3;
    orient(6, 2) = nz3;

    orient(7, 0) = nx4;
    orient(7, 1) = ny4;
    orient(7, 2) = nz4;

    orient(8, 0) = nx4;
    orient(8, 1) = ny4;
    orient(8, 2) = nz4;

    orient(9, 0) = nx4;
    orient(9, 1) = ny4;
    orient(9, 2) = -nz4;

    GeneralPatch c(vec1, numb, params, orient);

    cout << "created patch" << endl;

    //int n2 = 100;
    //double packing_fraction = 0.01;

    double l = cbrt(pi * (double)m2 / (6. * packing_fraction));


    cout << l << endl;

    Condensate A(l, n);

    cout << "created condensate" << endl;
    //TwoTetrahedral c(10.0, 1.4, pi / 4., 0.0, 1., pi / 6., 0.0, 1., pi / 6., 1000, 1000);

    // string filp = "/home/dino/Desktop/Chemistry/SimulationResults/ChemicalOscillator/sim-20-12-14-19:43:58/pos_beta=1_i=0455.csv";
    // string filo = "/home/dino/Desktop/Chemistry/SimulationResults/ChemicalOscillator/sim-20-12-14-19:43:58/orientation_beta=1_i=0455.csv";

    // string filp = "/home/dino/Documents/Condensate/TernaryFluid2/pos_beta=1_i=02097.csv";
    // string filo = "/home/dino/Documents/Condensate/TernaryFluid2/orientation_beta=1_i=02097.csv";

    // double T;
    // bool err1;
    // bool err2;
    // matrix<double> temppos = importcsv(filp, T, err1);
    // matrix<double> tempori = importcsv(filo, T, err1);

    // matrix<double> newpos(2000,3);
    // matrix<double> newori(2000,3);

    // A.obj->setdat(temppos);
    // A.obj->setorientation(tempori);

    // cout << b.tripr111 << endl;
    // cout << b.doubr11 << endl;
    // cout << b.doubr22 << endl;
    // cout << b2.on_rate << endl;
    // cout << b2.off_rate << endl;
    // pausel();

    //TetrahedralPatch c2(10.0,1.4,0.927);

    A.setBindingModel(b);

    //cout << "set up 1" << endl;
    A.setpots(c);

    //int a = system("python3 /home/dino/Documents/Condensate/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

    //int a = system("python3 /home/dino/Desktop/tylercollab/Repo/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

    //A.run_singlebond(10000, 1000);
    // for(int i = 0 ; i < 6 ; i++) {
    
    // for(int j = 0 ; j < (c.i1)[i][0] ; j++ ) {
    // cout << (c.i1)[i][j+1] <<" ";
    
    // }
    // cout << endl;
    // }

  

    A.setviscosity(0.1);

    double beta =1.0;

    A.obj->setkT(1. / beta);

    stringstream ss;
    ss << unbinding_rebar;

    string base = "_rebar=";
    base += ss.str();

    double T;
    int TT;
    bool vv1, vv2, vv3;

    matrix<double> postemp = importcsv("/u/home/d/dinoo/Chemistry/Code/Basic/InitialConditions/bigstartp.csv", T, vv1);
    matrix<double> orienttemp = importcsv("/u/home/d/dinoo/Chemistry/Code/Basic/InitialConditions/bigstarto.csv", T, vv2);
    matrix<int> bindtemp = importcsv("/u/home/d/dinoo/Chemistry/Code/Basic/InitialConditions/bigstartb.csv", TT, vv3);

    // matrix<double> postemp = importcsv("/home/dino/Desktop/tylercollab/Repo/Code/Basic/InitialConditions/bigstartp.csv", T, vv1);
    // matrix<double> orienttemp = importcsv("/home/dino/Desktop/tylercollab/Repo/Code/Basic/InitialConditions/bigstarto.csv", T, vv2);
    // matrix<int> bindtemp = importcsv("/home/dino/Desktop/tylercollab/Repo/Code/Basic/InitialConditions/bigstartb.csv", TT, vv3);

    A.obj->setdat(postemp);
    A.obj->setorientation(orienttemp);
    BinaryBindStore bbs2;
    vector1<bool> iss(bindtemp.getncols());
    vector1<int> ist(bindtemp.getncols());
    for (int i = 0; i < bindtemp.getncols(); i++)
    {
        iss[i] = (bool)bindtemp(0, i);
        ist[i] = bindtemp(1, i);
    }
    bbs2.isbound = iss;
    bbs2.boundto = ist;
    //Do processing to make sure everything is fine here

    A.run_singlebond_different_sizes_continue(runtime, 1000, m2, 0, bbs2, base);

    //A.run_singlebond_different_sizes(100000, 10,m2, base);

    // A.run_singlebond_different_sizes(10000000, 1000, 6000, base);

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