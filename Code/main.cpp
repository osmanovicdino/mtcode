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

    

    int runtime;
    double packing_fraction;
    double size_of_part;
    int num_arms;
    double angs;
    double energ;

    if (argc == 7)
    {
        runtime = atof(argv[1]);
        packing_fraction = atof(argv[2]);
        size_of_part = atof(argv[3]);
        energ = atof(argv[4]);
        angs = atof(argv[5]);
        num_arms = atof(argv[6]);
    }
    else
    {
        error("incorrect arg number");
    }

    int nt = 10000;
    double l = cbrt(pi * CUB(size_of_part)*(double)nt / (6. * packing_fraction));

    Condensate A(l, nt);

    A.size_mol = size_of_part;

    vector1<int> vec1(1);
    vec1[0] = num_arms;

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

    matrix<double> params2(tot*tot, 3);

    for (int i = 0; i < tot * tot; i++)
    {
        params2(i, 0) = energ;
        params2(i, 1) = 1.4*size_of_part;
        params2(i, 2) = angs;
    }


    matrix<double> orient(num_arms, 3);

    cout << orient << endl;

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

    double nx1_2 = sqrt(8. / 9.);
    double ny1_2 = 0.;
    double nz1_2 = -1. / 3.;

    double nx2_2 = -sqrt(2. / 9.);
    double ny2_2 = sqrt(2. / 3.);
    double nz2_2 = -1. / 3.;

    double nx3_2 = -sqrt(2. / 9.);
    double ny3_2 = -sqrt(2. / 3.);
    double nz3_2 = -1. / 3.;

    double nx4_2 = 0;
    double ny4_2 = 0;
    double nz4_2 = 1.;

    matrix<double> asd2(4, 3);

    asd2(0, 0) = nx1_2;
    asd2(0, 1) = ny1_2;
    asd2(0, 2) = nz1_2;

    asd2(1, 0) = nx2_2;
    asd2(1, 1) = ny2_2;
    asd2(1, 2) = nz2_2;

    asd2(2, 0) = nx3_2;
    asd2(2, 1) = ny3_2;
    asd2(2, 2) = nz3_2;

    asd2(3, 0) = nx4_2;
    asd2(3, 1) = ny4_2;
    asd2(3, 2) = nz4_2;

    int iter2 = 0;
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < num_arms; j++)
        {
            if(num_arms == 3) {
            orient(iter2, 0) = asd(j, 0);
            orient(iter2, 1) = asd(j, 1);
            orient(iter2, 2) = asd(j, 2);
            iter2++;
            }
            else if(num_arms == 4) {
            orient(iter2, 0) = asd(j, 0);
            orient(iter2, 1) = asd(j, 1);
            orient(iter2, 2) = asd(j, 2);
            iter2++;
            }
        }
    }


    GeneralPatch c4(vec1, numb, params2, orient);

    BindingModelSingle b2(0.99999, 0.0001);




    // vector1<bool> pb(3, false);
    // cube geo(l, pb, 3);
    //A.obj->setgeometry(geo);

    A.setBindingModel(b2);

    A.setpots(c4);

    

    //double beta = 1.;
    double beta = 1.0;

    A.obj->setkT(1. / beta);

    string base = "den=";
    stringstream dd;
    dd << packing_fraction;
    base += dd.str();

    base += "_size_of_part=";
    stringstream ss;
    ss << size_of_part;
    base += ss.str();

    stringstream ss1;
    ss1 << energ;
    base += "_energy=";
    base += ss1.str();

    base += "_angular=";
    stringstream ss2;
    ss2 << angs;
    base += ss2.str();

    base += "_arms=";
    stringstream ii4;
    ii4 << num_arms;
    base += ii4.str();

    A.setviscosity(0.1*size_of_part);

  
    A.run_singlebond(runtime, 1000, base);
    
    return 0;
}