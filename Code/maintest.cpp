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
    int n = 400;

    Condensate A(27.5, n);

    double angle = 0.5;
    
    int tot = 4 * 4 + 4 * 2 + 2 * 2;
    matrix<double> params(tot, 3);
    for (int i = 0; i < 16; i++)
    {
        params(i, 0) = 0.0;
        params(i, 1) = 1.4;
        params(i, 2) = angle;
    }
    for (int i = 16; i < 24; i++)
    {
        params(i, 0) = 100.0;
        params(i, 1) = 1.4;
        params(i, 2) = angle;
    }
    for (int i = 24; i < tot; i++)
    {
        params(i, 0) = 100.0;
        params(i, 1) = 1.4;
        params(i, 2) = angle;
    }

    TetrahedralWithBivalent c2(params, 50,n);
    A.setpots(c2);

    int runtime = 1000000;
    int every = 1000;
    A.setviscosity(10.0);
    // A.run_singlebond(runtime, every);

    A.run(10000000, 1000);

    return 0;
}