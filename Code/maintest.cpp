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
#include "Nanotube/Nanotube.h"
#include "DataStructures/importbmp.h"
#include "DataStructures/interpolation.h"
// #include "NCGasR.h"
// #include "Microtubule.h"

//#include "MDGPU.cu"

using namespace std;

void done() {
    cout << "done" << endl;
}

int main(int argc, char **argv)
{

    srand(time(NULL));

    int N = 1000;
    double pd = 0.01;
    double d = cbrt(N/pd);

    NanotubeAssembly A(d,N);

    // BivalentPatch p2(100.0,1.3,pi/3.);
    // A.setpots(p2);
    //A.setkT(0.1);


    A.run(1000000,1000);
/* 

    vector1<int> vec1(4);
    vec1[0]=4;
    vec1[1]=3;
    vec1[2]=2;
    vec1[3]=1;

    vector1<int> numb(4);

    numb[0] = 2000;
    numb[1] = 3000;
    numb[2] =  3500;
    numb[3] =  7500;

    int tot = 0;
    for(int i = 0 ; i < 4 ; i++) {
        for(int j = i ; j < 4 ; j++) {
            tot += vec1[i]*vec1[j];
        }
    }
    matrix<double> params(tot,3);

    for(int i = 0  ; i < tot ; i++) {
        params(i,0) = 10.0;
        params(i,1) = 1.4;
        params(i,2) = 0.927;
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

    matrix<double> asd(4,3);

    asd(0,0) = nx1;
    asd(0,1) = ny1;
    asd(0,2) = nz1;

    asd(1, 0) = nx2;
    asd(1, 1) = ny2;
    asd(1, 2) = nz2;

    asd(2, 0) = nx3;
    asd(2, 1) = ny3;
    asd(2, 2) = nz3;

    asd(3, 0) = nx4;
    asd(3, 1) = ny4;
    asd(3, 2) = nz4;

    int iter = 0;
    for(int i = 0  ; i < 4 ; i++) {
        for(int j = 0  ; j < vec1[i] ; j++) {
        orient(iter, 0) = asd(j, 0);
        orient(iter, 1) = asd(j, 1);
        orient(iter, 2) = asd(j, 2);
        iter++;
        }
    }

    GeneralPatch a(vec1,numb, params,orient);



    cout << a.no_patches_per_type << endl;
    cout << a.total_patches_per_type << endl;    
    cout << a.num_per_type << endl;

    int i,j;
    a.which_particle(11500,15786,i,j);
    cout << i << " " << j << endl;

    cout << a.return_type(i) << " " << a.return_type(j) << endl;

    int wpi,wpj;
    int potn = a.pot_starters(2,3);
    cout << potn << endl;
    a.which_patch(i,j,potn,wpi,wpj);

    cout << wpi << " " << wpj << endl;
    
    cout << a.which_potential(1000,3100,wpi,wpj) << endl;

    
    pausel();
     */


    // BMP a("./Basic/InitialConditions/TestImage.bmp");

    // Bilinear3 a(20., 20., "./Basic/InitialConditions/TestImage.bmp");


    // a.gamma_white = 2.;
    // a.gamma_black = 1.;

    // cout << a.data << endl;

    // for(int i = 0 ; i < 20 ; i++) {
    //     for(int j = 0 ; j < 20 ; j++) {
    //         cout << a(double(i),double(j)) << ", ";
    //     }
    //     cout << endl;
    // }

   // a(10,10);



    // cout << a.data.size() << endl;

    // matrix<int> b(a.bmp_info_header.width,a.bmp_info_header.height);

    // int iter = 0;
    // for(int i = 0  ; i < a.bmp_info_header.height ; i++) {
    //     for(int j = 0  ; j < a.bmp_info_header.width ; j++) {
    //         b(i,j) = a.data[iter];
    //         iter++;
    //     }
    // }







    // for(int i = 0 ; i < 28000 ; i++) {
    //     asd[i] = rand() % 28000;
    // }

    // for(int j = 0 ; j < 1000 ; j++) {
    //     sort(asd.begin(),asd.end());
    // }

}