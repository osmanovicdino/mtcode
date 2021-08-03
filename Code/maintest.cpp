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

#include "fftw3.h"

using namespace std;



int main(int argc, char **argv)
{

    int ind1[500] = { 2,
                      4,
                      5,
                      6,
                      6,
                      7,
                      7,
                      9,
                      9,
                      15,
                      16,
                      17,
                      18,
                      19,
                      20,
                      22,
                      23,
                      27,
                      28,
                      29,
                      30,
                      31,
                      31,
                      31,
                      33,
                      34,
                      37,
                      37,
                      37,
                      38,
                      39,
                      40,
                      40,
                      41,
                      41,
                      42,
                      42,
                      44,
                      45,
                      46,
                      46,
                      47,
                      48,
                      49,
                      51,
                      52,
                      53,
                      54,
                      55,
                      56,
                      56,
                      57,
                      58,
                      62,
                      64,
                      65,
                      66,
                      69,
                      69,
                      70,
                      70,
                      70,
                      73,
                      74,
                      74,
                      75,
                      76,
                      77,
                      80,
                      80,
                      80,
                      80,
                      81,
                      82,
                      84,
                      86,
                      86,
                      88,
                      88,
                      88,
                      89,
                      91,
                      93,
                      94,
                      94,
                      94,
                      95,
                      95,
                      95,
                      96,
                      96,
                      97,
                      97,
                      99,
                      102,
                      102,
                      104,
                      105,
                      107,
                      107,
                      107,
                      107,
                      108,
                      109,
                      111,
                      112,
                      113,
                      115,
                      119,
                      120,
                      123,
                      124,
                      124,
                      124,
                      125,
                      125,
                      130,
                      133,
                      133,
                      134,
                      137,
                      137,
                      138,
                      138,
                      139,
                      139,
                      140,
                      141,
                      143,
                      143,
                      144,
                      144,
                      145,
                      145,
                      146,
                      148,
                      150,
                      152,
                      155,
                      155,
                      157,
                      157,
                      158,
                      159,
                      160,
                      161,
                      162,
                      165,
                      165,
                      166,
                      167,
                      168,
                      169,
                      169,
                      169,
                      170,
                      171,
                      171,
                      172,
                      173,
                      173,
                      174,
                      174,
                      176,
                      178,
                      178,
                      183,
                      184,
                      185,
                      185,
                      186,
                      187,
                      187,
                      188,
                      191,
                      193,
                      195,
                      195,
                      195,
                      196,
                      197,
                      197,
                      198,
                      199,
                      200,
                      201,
                      202,
                      202,
                      204,
                      208,
                      209,
                      211,
                      213,
                      213,
                      213,
                      217,
                      222,
                      222,
                      223,
                      224,
                      225,
                      226,
                      226,
                      228,
                      229,
                      230,
                      233,
                      234,
                      235,
                      236,
                      237,
                      241,
                      242,
                      242,
                      242,
                      245,
                      247,
                      248,
                      249,
                      251,
                      251,
                      251,
                      253,
                      257,
                      258,
                      261,
                      262,
                      264,
                      265,
                      266,
                      270,
                      271,
                      271,
                      271,
                      273,
                      274,
                      275,
                      276,
                      277,
                      278,
                      279,
                      283,
                      284,
                      288,
                      290,
                      293,
                      295,
                      295,
                      298,
                      299,
                      299,
                      299,
                      300,
                      301,
                      307,
                      308,
                      309,
                      310,
                      310,
                      311,
                      312,
                      312,
                      312,
                      313,
                      314,
                      314,
                      319,
                      321,
                      323,
                      324,
                      328,
                      328,
                      330,
                      332,
                      336,
                      337,
                      337,
                      341,
                      343,
                      344,
                      345,
                      345,
                      346,
                      349,
                      353,
                      356,
                      359,
                      360,
                      362,
                      364,
                      364,
                      364,
                      371,
                      371,
                      372,
                      377,
                      378,
                      381,
                      381,
                      382,
                      384,
                      387,
                      387,
                      389,
                      389,
                      391,
                      392,
                      392,
                      395,
                      397,
                      397,
                      399,
                      402,
                      402,
                      406,
                      408,
                      409,
                      409,
                      412,
                      412,
                      412,
                      414,
                      418,
                      419,
                      420,
                      424,
                      426,
                      426,
                      431,
                      434,
                      436,
                      436,
                      438,
                      440,
                      442,
                      442,
                      443,
                      445,
                      446,
                      448,
                      449,
                      453,
                      454,
                      455,
                      456,
                      463,
                      467,
                      469,
                      470,
                      470,
                      470,
                      470,
                      471,
                      479,
                      482,
                      483,
                      485,
                      490,
                      493,
                      498,
                      500,
                      503,
                      505,
                      510,
                      511,
                      512,
                      512,
                      512,
                      512,
                      513,
                      514,
                      514,
                      516,
                      523,
                      524,
                      525,
                      525,
                      528,
                      529,
                      530,
                      533,
                      537,
                      540,
                      541,
                      542,
                      553,
                      555,
                      556,
                      557,
                      558,
                      558,
                      562,
                      563,
                      565,
                      566,
                      567,
                      567,
                      567,
                      568,
                      568,
                      575,
                      576,
                      580,
                      581,
                      582,
                      590,
                      591,
                      591,
                      595,
                      595,
                      596,
                      597,
                      600,
                      601,
                      609,
                      609,
                      611,
                      611,
                      614,
                      617,
                      618,
                      621,
                      621,
                      630,
                      639,
                      640,
                      641,
                      642,
                      645,
                      646,
                      649,
                      653,
                      654,
                      657,
                      661,
                      663,
                      663,
                      666,
                      667,
                      669,
                      669,
                      671,
                      671,
                      671,
                      671,
                      676,
                      682,
                      683,
                      684,
                      685,
                      692,
                      693,
                      698,
                      700,
                      702,
                      706,
                      709,
                      710,
                      710,
                      711,
                      714,
                      720,
                      721,
                      729,
                      731,
                      734,
                      735,
                      739,
                      744,
                      755,
                      761,
                      765,
                      769,
                      774,
                      777,
                      795,
                      796,
                      798,
                      805,
                      813,
                      819,
                      820,
                      820,
                      823,
                      826,
                      828,
                      829,
                      833,
                      839,
                      853,
                      877,
                      892,
                      902,
                      929,
                      930,
                      932,
                      943,
                      961,
                      962,
                      990 };

    int ind2[500] = { 789,
                      653,
                      905,
                      912,
                      18,
                      73,
                      465,
                      906,
                      73,
                      848,
                      25,
                      899,
                      275,
                      489,
                      631,
                      351,
                      699,
                      758,
                      144,
                      990,
                      218,
                      719,
                      937,
                      998,
                      734,
                      856,
                      928,
                      265,
                      307,
                      62,
                      702,
                      775,
                      210,
                      376,
                      598,
                      270,
                      271,
                      816,
                      118,
                      874,
                      968,
                      62,
                      208,
                      258,
                      200,
                      328,
                      189,
                      558,
                      980,
                      920,
                      665,
                      957,
                      441,
                      880,
                      532,
                      283,
                      102,
                      689,
                      526,
                      199,
                      249,
                      741,
                      575,
                      637,
                      103,
                      687,
                      98,
                      692,
                      285,
                      900,
                      457,
                      447,
                      586,
                      951,
                      567,
                      984,
                      217,
                      314,
                      995,
                      98,
                      208,
                      678,
                      899,
                      636,
                      640,
                      705,
                      221,
                      288,
                      585,
                      550,
                      589,
                      709,
                      231,
                      503,
                      654,
                      659,
                      130,
                      312,
                      987,
                      492,
                      666,
                      683,
                      711,
                      962,
                      463,
                      392,
                      671,
                      512,
                      809,
                      340,
                      635,
                      319,
                      786,
                      874,
                      134,
                      423,
                      946,
                      196,
                      977,
                      893,
                      268,
                      636,
                      379,
                      427,
                      786,
                      747,
                      197,
                      463,
                      502,
                      583,
                      832,
                      184,
                      963,
                      534,
                      721,
                      845,
                      957,
                      273,
                      332,
                      992,
                      551,
                      715,
                      859,
                      762,
                      418,
                      492,
                      785,
                      465,
                      489,
                      503,
                      695,
                      702,
                      481,
                      649,
                      593,
                      224,
                      430,
                      378,
                      671,
                      954,
                      978,
                      172,
                      247,
                      838,
                      295,
                      937,
                      971,
                      494,
                      660,
                      818,
                      820,
                      583,
                      609,
                      692,
                      845,
                      598,
                      855,
                      652,
                      701,
                      957,
                      943,
                      356,
                      622,
                      675,
                      703,
                      299,
                      339,
                      918,
                      223,
                      504,
                      466,
                      619,
                      935,
                      789,
                      636,
                      630,
                      720,
                      468,
                      719,
                      283,
                      471,
                      889,
                      655,
                      883,
                      982,
                      747,
                      832,
                      594,
                      993,
                      384,
                      239,
                      383,
                      805,
                      636,
                      798,
                      838,
                      431,
                      715,
                      566,
                      781,
                      935,
                      568,
                      670,
                      811,
                      711,
                      308,
                      996,
                      833,
                      973,
                      965,
                      343,
                      918,
                      652,
                      656,
                      492,
                      664,
                      706,
                      814,
                      982,
                      274,
                      509,
                      881,
                      896,
                      309,
                      922,
                      636,
                      638,
                      988,
                      820,
                      976,
                      380,
                      472,
                      318,
                      352,
                      361,
                      478,
                      679,
                      770,
                      509,
                      970,
                      437,
                      938,
                      606,
                      647,
                      964,
                      972,
                      331,
                      760,
                      447,
                      628,
                      940,
                      360,
                      635,
                      951,
                      868,
                      701,
                      346,
                      900,
                      326,
                      844,
                      521,
                      649,
                      499,
                      656,
                      684,
                      592,
                      717,
                      733,
                      816,
                      362,
                      748,
                      782,
                      860,
                      597,
                      536,
                      337,
                      319,
                      404,
                      819,
                      639,
                      868,
                      527,
                      714,
                      967,
                      456,
                      340,
                      588,
                      379,
                      912,
                      419,
                      468,
                      746,
                      360,
                      861,
                      869,
                      630,
                      679,
                      901,
                      492,
                      746,
                      986,
                      480,
                      989,
                      783,
                      838,
                      994,
                      475,
                      469,
                      883,
                      381,
                      544,
                      582,
                      628,
                      696,
                      541,
                      686,
                      488,
                      675,
                      885,
                      930,
                      695,
                      702,
                      833,
                      785,
                      835,
                      989,
                      906,
                      701,
                      861,
                      557,
                      708,
                      551,
                      567,
                      772,
                      544,
                      629,
                      667,
                      796,
                      797,
                      635,
                      912,
                      581,
                      704,
                      626,
                      966,
                      758,
                      727,
                      953,
                      695,
                      461,
                      743,
                      781,
                      913,
                      629,
                      993,
                      519,
                      593,
                      595,
                      760,
                      644,
                      743,
                      928,
                      479,
                      487,
                      680,
                      962,
                      578,
                      620,
                      720,
                      557,
                      734,
                      969,
                      601,
                      733,
                      907,
                      721,
                      757,
                      515,
                      924,
                      826,
                      700,
                      776,
                      787,
                      959,
                      857,
                      567,
                      548,
                      834,
                      915,
                      692,
                      721,
                      796,
                      813,
                      997,
                      870,
                      574,
                      940,
                      728,
                      792,
                      657,
                      673,
                      846,
                      969,
                      734,
                      976,
                      945,
                      784,
                      909,
                      929,
                      716,
                      865,
                      706,
                      716,
                      872,
                      720,
                      677,
                      846,
                      953,
                      675,
                      808,
                      987,
                      901,
                      907,
                      925,
                      864,
                      704,
                      715,
                      954,
                      974,
                      744,
                      855,
                      921,
                      937,
                      938,
                      947,
                      855,
                      990,
                      750,
                      825,
                      716,
                      721,
                      904,
                      859,
                      983,
                      759,
                      781,
                      965,
                      866,
                      980,
                      882,
                      825,
                      904,
                      891,
                      905,
                      940,
                      839,
                      982,
                      873,
                      871,
                      930,
                      911,
                      982,
                      868,
                      831,
                      879,
                      903,
                      912,
                      954,
                      999,
                      976,
                      947,
                      992,
                      958,
                      965,
                      955,
                      964,
                      993,
                      975,
                      997,
                      982 };

    vector<mdpairwd> edgelist(500);
    for(int i = 0  ; i < 500 ; i++) {
        int wp1 =  ind1[i];
        int wp2 =  ind2[i];
        double en = 0.0;
        mdpairwd test(wp1, wp2, en); //now our score is the energy
        edgelist.push_back(test);
    }

    string sg = "a";
    vector1<int> indexes2(1000, sg);

    ConnectedComponentsParallel(edgelist,indexes2);

    cout << indexes2 << endl;
    
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