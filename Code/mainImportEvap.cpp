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

    srand(time(NULL));
    //signal(SIGSEGV, handler);

    double int1;
    double int2;
    double int3;
    double int4;
    int runtime;
    double energy_barrier;
    if (argc == 7)
    {
        runtime = atof(argv[1]);
        int1 = atof(argv[3]);
        int2 = atof(argv[4]);
        int3 = atof(argv[5]);
        int4 = atof(argv[6]);
        energy_barrier = atof(argv[9]);
    }
    else
    {
        //error("incorrect arg number");
        runtime = 1000001;
        int1 = 10.0;
        int2 = 10.0;
        int3 = 10.0;
        int4 = 50.0;
        energy_barrier = 0.01;
    }

    cout << " " << int1 << " " << int2 << "  " << int3 << endl;

    // for(int i = 0 ; i < 28 ; i++) {
    //     params(i, 0) =  10.0; //strength
    //     params(i, 1) =  1.4; //distance
    //     params(i, 2) = 0.927;
    // }
    int m1 = 1000;
    int m2 = m1 + 0;
    int n = m2 + 2000;

    // int m1  = 2;
    // int m2 = 6;
    // int n = 10;

    //SortingFunctionUniform my_sorter;

    SortingFunctionNonUniform my_sorter;
    my_sorter.div1 = (m1 * 4);
    my_sorter.div2 = (m1 * 4 + 4 * (m2 - m1));
    //BindingModelTernary b(m1 * 4, m1*4 + 4 * (m2-m1));

    BindingModelTernary<SortingFunctionNonUniform> b(my_sorter);
    //BindingModelTernary b(my_sorter);

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

    double unbinding_rebar = -20.0;
    double unbinding_rebar2 = 2.0;

    b.setup_energy_barrier(0.99999, 0.99999, 0.99999, 0.99999, energy_barrier, 0.99999,
                           0.99999, 0.99999, 0.99999, 0.99999, 0.99999, 0.99999,
                           0.0,
                           unbinding_rebar2,
                           0.0,
                           0.0,
                           0.0,
                           unbinding_rebar,
                           0.0,
                           0.0,
                           0.0);

    bool c1;

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

    double size1 = 2.0;
    double size2 = 1.0;

    double sizemix = (size1 + size2) / 2.;
    // pausel();

    vector1<int> vec1(3);
    vec1[0] = 4;
    vec1[1] = 4;
    vec1[2] = 3;

    vector1<int> numb(3);

    numb[0] = m1;
    numb[1] = m2;
    numb[2] = n;

    int tot = 4 * 4 + 4 * 4 + 4 * 3 + 4 * 4 + 4 * 3 + 3 * 3;
    matrix<double> params(tot, 3);

    int iter = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            params(iter, 0) = int1;
            params(iter, 1) = 1.4 * size1;
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
                params(iter, 1) = 1.4 * size1;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1.4 * size1;
                params(iter, 2) = 0.927;
            }

            iter++;
        }
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (j == 0)
            {
                params(iter, 0) = int2;
                params(iter, 1) = 1.4 * sizemix;
                params(iter, 2) = 0.927;
            }
            else if (j == 1)
            {
                params(iter, 0) = int3;
                params(iter, 1) = 1.4 * sizemix;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = int4;
                params(iter, 1) = 1.4 * sizemix;
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
                params(iter, 1) = 1 * size1;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1. * size1;
                params(iter, 2) = 0.927;
            }

            iter++;
        }
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
        {

            if (j == 2)
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 2. * sizemix; //destroy the gas phase
                params(iter, 2) = 0.927;
                iter++;
            }
            else
            {

                params(iter, 0) = 0.0;
                params(iter, 1) = 1.4 * sizemix;
                params(iter, 2) = 0.927;
                iter++;
            }
        }
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i == 1 && j == 1)
            {
                params(iter, 0) = int3;
                params(iter, 1) = 1.4 * size2;
                params(iter, 2) = 0.927;
            }
            else
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = 1.4 * size2;
                params(iter, 2) = 0.927;
            }
            iter++;
        }
    }

    matrix<double> orient(11, 3);

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

    double nx5 = sqrt(3.) / 2.;
    double ny5 = -1. / 2.;
    double nz5 = 0.;

    double nx6 = 0.;
    double ny6 = 1.;
    double nz6 = 0.;

    double nx7 = -sqrt(3.) / 2.;
    double ny7 = -1. / 2.;
    double nz7 = 0.;

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

    orient(8, 0) = nx5;
    orient(8, 1) = ny5;
    orient(8, 2) = nz5;

    orient(9, 0) = nx6;
    orient(9, 1) = ny6;
    orient(9, 2) = nz6;

    orient(10, 0) = nx7;
    orient(10, 1) = nz7;
    orient(10, 2) = ny7;

    GeneralPatch c(vec1, numb, params, orient);

    cout << "created patch" << endl;

    //int n2 = 100;
    //double packing_fraction = 0.01;

    double l = cbrt(pi * CUB(size1) * (double)m1 / (6. * 0.02));

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

    double beta = 1.0;

    A.obj->setkT(1. / beta);

    stringstream ss;
    ss << unbinding_rebar;

    string base = "_rebar=";
    base += ss.str();

    // cout << "done" << endl;
    //Do processing to make sure everything is fine here

    //pausel();
    //cout << m2 << endl;
    A.obj->setdt(0.005);

    //A.run_singlebond(runtime, 1000, base);
    // vector<string> orientfiles;
    // vector<string> posfiles;
    // vector<string> bindfiles;

    // return_csv_in_current_dir("orient", orientfiles);
    // return_csv_in_current_dir("pos", posfiles);
    // return_csv_in_current_dir("bind", bindfiles);

    // int s1 = orientfiles.size();
    // int s2 = posfiles.size();
    // int s3 = bindfiles.size();

    // if ((s1 == 0) || (s2 == 0) || (s3 == 0)) {
    //     error("no files found to import");
    // }

    // if( (s1 != s2 ) || (s2 != s3 ) || (s1 != s3 ) ) {
    //     error("different sizes of import");
    // }

    double T;
    int TT;
    bool vv1,vv2,vv3;
    // matrix<double> postemp = importcsv("/home/dino/Desktop/tylercollab/Repo/Code/Basic/InitialConditions/posi.csv", T, vv1);
    // matrix<double> orienttemp = importcsv("/home/dino/Desktop/tylercollab/Repo/Code/Basic/InitialConditions/poso.csv", T, vv2);
    // matrix<int> bindtemp = importcsv("/home/dino/Desktop/tylercollab/Repo/Code/Basic/InitialConditions/posb.csv", TT, vv3);


    string dir = "/u/home/d/dinoo/Chemistry/Code/Basic/InitialConditions/";
    string dir2 = "/home/dino/Desktop/tylercollab/Repo/Code/Basic/InitialConditions/";

    matrix<double> postemp = importcsv(dir2+"posi.csv", T, vv1);
    matrix<double> orienttemp = importcsv(dir2+"poso.csv", T, vv2);
    matrix<int> bindtemp = importcsv(dir2+"posb.csv", TT, vv3);


    A.obj->setdat(postemp);
    A.obj->setorientation(orienttemp);
    BinaryBindStore bbs2;
    vector1<bool> iss(bindtemp.getncols());
    vector1<int> ist(bindtemp.getncols());
    for(int i  = 0 ; i < bindtemp.getncols() ; i++ ) {
        iss[i] = (bool)bindtemp(0,i);
        ist[i] =  bindtemp(1,i);
    }
    bbs2.isbound =  iss;
    bbs2.boundto = ist;
    //Do processing to make sure everything is fine here

    int every = 1000;

    A.run_singlebond_different_sizes_continue(runtime, every, m2, size1, size2, 0, bbs2, base);

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