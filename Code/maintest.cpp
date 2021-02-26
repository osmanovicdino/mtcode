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

    // BMP a("./Basic/InitialConditions/TestImage.bmp");

    Bilinear3 a(20., 20., "./Basic/InitialConditions/TestImage.bmp");


    a.gamma_white = 2.;
    a.gamma_black = 1.;

    cout << a.data << endl;

    for(int i = 0 ; i < 20 ; i++) {
        for(int j = 0 ; j < 20 ; j++) {
            cout << a(double(i),double(j)) << ", ";
        }
        cout << endl;
    }

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