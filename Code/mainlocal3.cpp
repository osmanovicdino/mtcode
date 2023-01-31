#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
// #include <limits>
// #include <cmath>
// #include <complex>
#include <sstream>
#include <string>
#include <iomanip>
// #include <sys/ioctl.h>
// #include <fcntl.h>
// #include <time.h>
// #include <sys/time.h>
#include <sys/stat.h>
#include <random>
#include <algorithm>
#include <parallel/algorithm>
#include <chrono>

//includes needed for backtracing

#include <unistd.h>
#include <execinfo.h>
#include <signal.h>

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

using namespace std;
using namespace std::chrono;

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



int main(int argc, char **argv)
{

    srand(time(NULL));
    signal(SIGSEGV, handler);



    




    
    double packing_fraction;
    double int1;
    double int2;
    double int3;
    double int4;
    double int5;
    int m3;
    int m4;

    int runtime;
    double energy_barrier;
    double anti_en;
    double inv_en;
    if (argc == 7)
    {
        runtime = atof(argv[1]);
        packing_fraction = atof(argv[2]);
        int1 = atof(argv[3]);
        int2 = atof(argv[4]);
        m3 = atof(argv[5]);
        m4 = atof(argv[6]);
    }
    else
    {
        //error("incorrect arg number");
        //error("incorrect argument number");
        runtime = 1000001;
        packing_fraction = 0.005;
        int1 = 15.0;
        int2 = 0.0;
        m3 = 500;
        m4 = 500;

    }

    //cout << packing_fraction << " " << int1 << " " << int2 << "  " << int3 << endl;

    // for(int i = 0 ; i < 28 ; i++) {
    //     params(i, 0) =  10.0; //strength
    //     params(i, 1) =  1.4; //distance
    //     params(i, 2) = 0.927;
    // }
    int m1 = m3;
    int m2 = m3+m4;


    // int m1  = 2;
    // int m2 = 6;
    // int n = 10;

    //SortingFunctionUniform my_sorter;

    // SortingFunctionNonUniform my_sorter;
    // my_sorter.div1 = (m1 * 4);
    // my_sorter.div2 = (m1 * 4 + 4 * (m2 - m1));
    // my_sorter.np = 4;

    //BindingModelTernary b(m1 * 4, m1*4 + 4 * (m2-m1));




    // BindingModelTernary<SortingFunctionNonUniform> b(my_sorter);
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

    BindingModelSingle b(0.999,0.001); 


    // double sizemix = (size1+size2)/2.;
    // pausel();

    double size = 1.0;
    int p1 = 3;
    int p2 = 3;
    vector1<int> vec1(2);
    vec1[0] = p1;
    vec1[1] = p2;

    vector1<int> numb(2);

    numb[0] = m1;
    numb[1] = m2;


    int tot   = p1 * p1 + p1 * p2 + p2 * p2;
    matrix<double> params(tot, 3);

    double range  = 1.4;
    double ang =  0.927;
    int iter = 0;
    for (int i = 0; i < p1; i++) 
    {
        for (int j = 0; j < p1; j++)
        {
            params(iter, 0) = int1;
            params(iter, 1) = range * size;
            params(iter, 2) = ang;
            iter++;
        }
    }

    for (int i = 0; i < p1; i++) 
    {
        for (int j = 0; j < p2; j++)
        {

            if (j == p2-1) // top one
            {
                params(iter, 0) = int2;
                params(iter, 1) = range * size;
                params(iter, 2) = ang; // slightly smaller aperture
            }

            else
            {
                params(iter, 0) = int1;
                params(iter, 1) = range * size;
                params(iter, 2) = ang; // slightly smaller aperture
            }
            iter++;
        }
    }

    for (int i = 0; i < p2; i++) 
    {
        for (int j = 0; j < p2; j++)
        {
            if (i==p1-1 || j == p2-1) //top one
            {
                params(iter, 0) = int2;
                params(iter, 1) = range * size;
                params(iter, 2) = ang; //slightly smaller aperture
            }

            else
            {
                params(iter, 0) = int1;
                params(iter, 1) = range * size;
                params(iter, 2) = ang; //slightly smaller aperture
            }
            iter++;
        }
    }


    matrix<double> orient(p1+p2, 3);

    // double nx1 = sqrt(8. / 9.);
    // double ny1 = 0.;
    // double nz1 = -1. / 3.;

    // double nx2 = -sqrt(2. / 9.);
    // double ny2 = sqrt(2. / 3.);
    // double nz2 = -1. / 3.;

    // double nx3 = -sqrt(2. / 9.);
    // double ny3 = -sqrt(2. / 3.);
    // double nz3 = -1. / 3.;



    double nx1 = 1.;
    double ny1 = 0.;
    double nz1 = 0.;

    double nx2 = -0.5;
    double ny2 = 0.5*sqrt(3.);
    double nz2 = 0.;

    double nx3 = -0.5 ;
    double ny3 = -0.5*sqrt(3.);
    double nz3 = 0.;

    // double nx8 =  0.5;
    // double ny8 = 0.5*sqrt(3);
    // double nz8 = 0.0;

    // double nx5 = 0.7100399393804211;
    // double ny5 = 0.7041614051370949;
    // double nz5 = 0.;

    // double nx6 = -0.9648416349034807;
    // double ny6 = 0.2628319226364604;
    // double nz6 = 0.;

    // double nx7 = 0.25480169552305926;
    // double ny7 = -0.9669933277735551;
    // double nz7 = 0.;

    // double nx8 = -0.2548016955230598;
    // double ny8 = 0.966993327773555;
    // double nz8 = 0.0;

    orient(0, 0) = nx1;
    orient(0, 1) = ny1;
    orient(0, 2) = nz1;

    orient(1, 0) = nx2;
    orient(1, 1) = ny2;
    orient(1, 2) = nz2;

    orient(2, 0) = nx3;
    orient(2, 1) = ny3;
    orient(2, 2) = nz3;

    orient(3, 0) = nx1;
    orient(3, 1) = ny1;
    orient(3, 2) = nz1;

    orient(4, 0) = nx2;
    orient(4, 1) = ny2;
    orient(4, 2) = nz2;

    orient(5, 0) = nx3;
    orient(5, 1) = ny3;
    orient(5, 2) = nz3;



    // orient(8, 0) = nx5;
    // orient(8, 1) = ny5;
    // orient(8, 2) = nz5;

    // orient(9, 0) = nx6;
    // orient(9, 1) = ny6;
    // orient(9, 2) = nz6;

    // orient(10, 0) = nx7;
    // orient(10, 1) = ny7;
    // orient(10, 2) = nz7;

    // orient(11, 0) = nx8;
    // orient(11, 1) = ny8;
    // orient(11, 2) = nz8;

    // orient(8, 0) = nx5;
    // orient(8, 1) = ny5;
    // orient(8, 2) = nz5;

    // orient(9, 0) = nx6;
    // orient(9, 1) = ny6;
    // orient(9, 2) = nz6;

    // orient(10, 0) = nx7;
    // orient(10, 1) = ny7;
    // orient(10, 2) = nz7;

    // orient(11, 0) = nx8;
    // orient(11, 1) = ny8;
    // orient(11, 2) = nz8;

    // orient(11, 0) = 0.;
    // orient(11, 1) = 0.;
    // orient(11, 2) = 1.;

    GeneralPatch c(vec1, numb, params, orient);

   // cout << "created patch" << endl;

    //TwoTetrahedral c2(10., 1.4, 0.927, 10., 1.4, 0.927, 10., 1.4, 0.927, 1000, 2000);
        //int n2 = 100;
        //double packing_fraction = 0.01;



    double l = cbrt(pi * CUB(size) * (double)(m1) / (6. * packing_fraction));

    //l = 20.34;

    Condensate A(l, m2 );


    A.setup_large_droplet(m3,m4,0,0.640658,  l);
    // outfunc(A.obj->getdat(),"res");

    //A.setup_tight_packing(1.8*size);

    
    //outfunc(A.obj->getdat(),"dat");
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

    string base = "den=";
    stringstream dd;
    dd << packing_fraction;
    base+= dd.str();
    
    
    base += "_int1=";
    stringstream ss;
    ss << int1;
    base += ss.str();

    base += "_int2=";
    stringstream ss2;
    ss2 << int2;
    base += ss2.str();

    // cout << "done" << endl;
    //Do processing to make sure everything is fine here
    base += "num_s1=";
    stringstream ss4;
    ss4 << m1;
    base += ss4.str();

    base += "num_s2=";
    stringstream ss5;
    ss5 << m2;
    base += ss5.str();

    // cout << "done" << endl;
    //Do processing to make sure everything is fine here

    //pausel();
    //cout << m2 << endl;
    A.obj->setdt(0.005);
    int every = 1000;

    //cout << "bout to run" << endl;

    // auto start = high_resolution_clock::now();



    A.run_singlebond(runtime, every, base);

    // auto stop = high_resolution_clock::now();

    // auto duration = duration_cast<microseconds>(stop - start);

    // cout << "Time taken by function: "
    //      << duration.count() << " microseconds" << endl;

    // BinaryBindStore bbs2;
    // int vak =  n*4;
    // vector1<bool> iss(vak);
    // vector1<int> ist(vak);
    // for (int i = 0; i < vak; i++)
    // {
    //     iss[i] = (bool)0;
    //     ist[i] = 0;
    // }
    // bbs2.isbound = iss;
    // bbs2.boundto = ist;

    // cout << "start" << endl;
    // A.run_singlebond_different_sizes_continue_thetalist(runtime, every, 1, size, size, 0, bbs2, base);
 

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