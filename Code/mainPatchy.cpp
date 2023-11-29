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


    double patchsize;
    
    if (argc == 2)
    {
        patchsize = atof(argv[1]);
    }





    
    int totN =2;
    double l = 5.;
    double int1 = 15.;
    int runtime = 10000000;
    int number_of_arms = 3;


    int m1 = totN;



    // int m1  = 2;
    // int m2 = 6;
    // int n = 10;

    //SortingFunctionUniform my_sorter;

    //BindingModelTernary b(m1 * 4, m1*4 + 4 * (m2-m1));
    double onrate=  0.99999;

    BindingModelSingle b(onrate, 1. - onrate);




    double size = 1.0; //size of the particle

    vector1<int> vec1(1);
    vec1[0] = number_of_arms;

    vector1<int> numb(1);

    numb[0] = m1;


    int tot = 3 * 3;
    matrix<double> params(tot, 3);

    double range  = 1.2; // range of the attraction
    // double patchsize = 0.927; //patch size in radians
    int iter = 0;
    for (int i = 0; i < number_of_arms; i++) //nanostar/nanostar interaction
    {
        for (int j = 0; j < number_of_arms; j++)
        {
            params(iter, 0) = int1;
            params(iter, 1) = range * size;
            params(iter, 2) = patchsize; //max angle of the patch
            iter++;
        }
    }




    matrix<double> orient(3, 3);


    // //DEFINE THE VECTORS TO THE PATCH ENDPOINTS
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
    double ny2 = sqrt(3. / 4.);
    double nz2 = 0.;

    double nx3 = -0.5;
    double ny3 = -sqrt(3. / 4.);
    double nz3 = 0.;

    orient(0, 0) = nx1;
    orient(0, 1) = ny1;
    orient(0, 2) = nz1;

    orient(1, 0) = nx2;
    orient(1, 1) = ny2;
    orient(1, 2) = nz2;

    orient(2, 0) = nx3;
    orient(2, 1) = ny3;
    orient(2, 2) = nz3;


    // orient(11, 0) = 0.;
    // orient(11, 1) = 0.;
    // orient(11, 2) = 1.;

    GeneralPatch c(vec1, numb, params, orient);

   // cout << "created patch" << endl;

    //TwoTetrahedral c2(10., 1.4, 0.927, 10., 1.4, 0.927, 10., 1.4, 0.927, 1000, 2000);
        //int n2 = 100;
        //double packing_fraction = 0.01;



    //l = 20.34;

    Condensate A(l, totN);


    
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

    A.setviscosity(10.);

    double beta = 1.0; //1/Temperature

    A.obj->setkT(1. / beta);

    string base = "l=";
    stringstream dd;
    dd << l;
    base+= dd.str();
    
    
    base += "_int1=";
    stringstream ss;
    ss << int1;
    base += ss.str();

    stringstream ss1;
    ss1 << patchsize;
    base += "_patchsize=";
    base += ss1.str();

    base += "_temperature=";
    stringstream ss2;
    ss2 << 1./beta;
    base += ss2.str();


    // cout << "done" << endl;
    //Do processing to make sure everything is fine here

    //pausel();
    //cout << m2 << endl;
    A.obj->setdt(0.005);
    int every = 1000;

    //cout << "bout to run" << endl;

    auto start = high_resolution_clock::now();



    A.run_singlebond(runtime, every, base);

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);

    cout << "Time taken by function: "
         << duration.count() << " microseconds" << endl;

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