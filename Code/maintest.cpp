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
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>

#include <unistd.h>
#include <execinfo.h>
#include <signal.h>


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
#include "MDBase/Analysis/AnalysisFunctions.h"
// #include "NCGasR.h"
// #include "Microtubule.h"

//#include "MDGPU.cu"

#include "fftw3.h"

using namespace std;



int main(int argc, char **argv)
{
    double l = 74.822;
    double T;
    bool err;
    //matrix<double> pos = importcsv("/home/dino/External/PhaseDiagramDesign6/den=0.005_i1=12._i2=12._i3=4._i4=60._m1=4000_m2=16000_rate=0.001/posden=0.005_int1=12_int2=12_int3=4_int4=60_br=0.001num_anti=4000num_inv=16000_i=02821.csv",T,err);
    matrix<double> pos = importcsv("/home/dino/External/PhaseDiagramDesign6/den=0.005_i1=12._i2=12._i3=4._i4=60._m1=16000_m2=16000_rate=0.001/posden=0.005_int1=12_int2=12_int3=4_int4=60_br=0.001num_anti=16000num_inv=16000_i=01865.csv", T, err);

    vector1<bool> pb(3, true);
    cube geo(l, pb, 3);

    // LangevinNVT *A;
    // A = new LangevinNVT(geo);

    double dis = 3.5;
    int num = floor(l / dis);
    int ccc;
    matrix<int> boxes = geo.generate_boxes_relationships(num, ccc);

    index_test * g = new index_test;
    vector1<int> indexes2 = distance_graph(pos, boxes,geo, 1.4, g);

    outfunc(indexes2, "myval4");

    // fftw_complex *in;
    // fftw_complex *out;

    // p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);



    /*
    int N = 6;
    double pd = 0.01;
    double d = cbrt(N/pd);

    NanotubeAssembly A(d,N);
    */
    
    // matrix<double> vf(2,3);
    // vf(0, 0) = d / 2.;
    // vf(0, 1) = d / 2.;
    // vf(0, 2) = d / 2.;

    // vf(1, 0) = d / 2. + 1.4;
    // vf(1, 1) = d / 2.;
    // vf(1, 2) = d / 2.;

    // matrix<double> ori2(2,9);


    // double theta =  0.0;
    // ori2(0, 0) = cos(theta);
    // ori2(0, 1) = 0.0;
    // ori2(0, 2) = sin(theta);
    // ori2(0, 3) = 0.0;
    // ori2(0, 4) = 1.0;
    // ori2(0, 5) = 0.0;
    // ori2(0, 6) = -sin(theta);
    // ori2(0, 7) = 0.0;
    // ori2(0, 8) = cos(theta);


    // double theta2 = 0.0;
    // ori2(1, 0) = -cos(theta2);
    // ori2(1, 1) = 0.0;
    // ori2(1, 2) = -sin(theta2);
    // ori2(1, 3) = 0.0;
    // ori2(1, 4) = 1.0;
    // ori2(1, 5) = 0.0;
    // ori2(1, 6) = sin(theta2);
    // ori2(1, 7) = 0.0;
    // ori2(1, 8) = -cos(theta2);

    // A.obj->setdat(vf);
    // A.obj->setorientation(ori2);
    // A.setkT(1.0);

    // A.run(1000,1);
    


    // BivalentPatch p2(100.0,1.3,pi/3.);
    // A.setpots(p2);
    //A.setkT(0.1);

    /*
    double beta = 1.;
    //for(double beta = 0.2 ; beta < 1.01 ; beta += 0.1) {
    double strp1[3] = {10.,20.,30.};
    double srang1[3] = {1.2, 1.3 , 1.4};
    double sang1[4] = {0.3,0.4,0.5,0.6};
    double sang2[4] = {0.1,0.2,0.3,0.4};

    double sdphi[4] = {0.1,0.2,0.3,0.4};
    double sdtheta[2] = {0.1,0.2};
    double sbaseangle[2] = {2 *pid / 2.5 , 2 *pid / 3.};
    




    for ( int i1 = 0 ; i1 < 3 ; i1++ )
    for ( int i2 = 0 ; i2 < 3;  i2++ )
    for (int i3 = 0 ; i3 < 4 ; i3++ )
    for (int i4 = 0 ; i4 < 3 ; i4++ )
    for (int i5 = 0 ; i5 < 3 ; i5++ )
    for (int i6 = 0 ; i6 < 4 ; i6++ )
    for (int i7 = 0  ; i7 < 4 ; i7++ )
    for (int i8 = 0  ; i8 < 2 ; i8++ )
    for (int i9 = 0  ; i9 < 2 ; i9++ ) {

    double str1 = strp1[i1];
    double rang1 = srang1[i2];
    double ang1 =  sang1[i3];

    double str2 = -strp1[i4];
    double rang2 = srang1[i5];
    double ang2 = sang2[i6];

    double dphi = sdphi[i7];
    double dtheta = sdtheta[i8];

    double baseangle = sbaseangle[i9];


    GeneralPatch c(CreateHexatic(N, str1, rang1, ang1, str2, rang2, ang2, dphi, dtheta, baseangle));
    stringstream ss;
    ss << i1;
    ss << i2;
    ss << i3;
    ss << i4;
    ss << i5;
    ss << i6;
    ss << i7;
    ss << i8;
    ss << i9;

    string base = "_identitystring=";
    base += ss.str();

    A.setkT(1. / beta);
    A.setpots(c);

    A.run(100000, 1000, base);
    }
    */

    //}
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







}