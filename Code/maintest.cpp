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

    matrix<double> *mom = new matrix<double>;
    matrix<double> *angmom = new matrix<double>;

    matrix<double> h(28000,3);

    matrix<double> F(28000,3,1.);

    matrix<double> T(28000,3,4.235457456746);

    *mom = h;
    *angmom = h;

    double dt  = 0.005;

    for(int j =0 ; j < 10000 ; j++) {
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < (*mom).getNsafe(); i++)
        {
            (*mom)(i, 0) += (dt / 2.) * F(i, 0);
            (*angmom)(i, 0) += (dt / 2.) * T(i, 0);
            (*mom)(i, 1) += (dt / 2.) * F(i, 1);
            (*angmom)(i, 1) += (dt / 2.) * T(i, 1);
            (*mom)(i, 2) += (dt / 2.) * F(i, 2);
            (*angmom)(i, 2) += (dt / 2.) * T(i, 2);
        }
    }

    done();


    // for(int i = 0 ; i < 28000 ; i++) {
    //     asd[i] = rand() % 28000;
    // }

    // for(int j = 0 ; j < 1000 ; j++) {
    //     sort(asd.begin(),asd.end());
    // }

}