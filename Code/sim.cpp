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
#include <unistd.h>
#include <execinfo.h>
#include <signal.h>
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
#include "Condensate/Nanostar.h"
// #include "NCGasR.h"
// #include "Microtubule.h"

//#include "MDGPU.cu"

using namespace std;

int main(int argc, char **argv)
{

    uint64_t microseconds_since_epoch = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // cout << microseconds_since_epoch << endl;
    // cout << time(NULL) << endl;
    int seed = microseconds_since_epoch % time(NULL);
    srand(seed);

    double dens;
    int num;
    if (argc == 3)
    {
        num = atof(argv[1]);
        dens = atof(argv[2]);
    }
    else
    {
        error("not enough arguments");
    }

    vector1<int> length_per_arm(4);
    length_per_arm[0] = 20;
    length_per_arm[1] = 10;
    length_per_arm[2] = 20;
    length_per_arm[3] = 10;

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

    matrix<double> ori(4,3);
    ori(0, 0) = nx1;
    ori(0, 1) = ny1;
    ori(0, 2) = nz1;

    ori(1, 0) = nx2;
    ori(1, 1) = ny2;
    ori(1, 2) = nz2;

    ori(2, 0) = nx3;
    ori(2, 1) = ny3;
    ori(2, 2) = nz3;

    ori(3, 0) = nx4;
    ori(3, 1) = ny4;
    ori(3, 2) = nz4;

    double ll =  cbrt((double)num/dens);

    Nanostar A(num,ll,length_per_arm,ori);


    A.DoAnMC();


    cout << "everything ?seems to work" << endl;

    //A.Passa_set_nanostar();

    A.run(1000000,1000);

    

    // string filename ="./";

    // int a = system("python3 /home/dino/Documents/Condensate/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

    // A.set_initial_state(filename);

    // for(int i = 0  ; i < A.bindpairs.size() ; i++) {
    // mdpair temp = A.bindpairs[i];
    // cout << temp.a << " " << temp.b << endl;
    // }

    // A.run(10000,100);
    return 0;

}