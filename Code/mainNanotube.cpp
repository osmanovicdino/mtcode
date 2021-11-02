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
#include "Nanotube/Nanotube.h"

// #include "NCGasR.h"
// #include "Microtubule.h"

//#include "MDGPU.cu"

using namespace std;

int main(int argc, char **argv)
{

    srand(time(NULL));
    int NM;

    if(argc == 2) {
        NM = atof(argv[1]);
    }

    //signal(SIGSEGV, handler);
    double radius = 25.;
    double monomers = NM;

    //monomers/(4/3piradius3)

    cout << "starting" << endl;
    NanotubeAssembly A(radius, monomers);

    double deltaG = 20.0;
    double angle = 0.6;
    BivalentPatch c2(deltaG, 1.4, angle);

    A.setpots(c2);
    A.setkT(1.0);

    cout << "done" << endl;
    ShellProperties B;
    int T;
    bool err1;
    matrix<int> pairs = importcsv("/home/dino/Documents/tylercollab/Repo/Code/Basic/InitialConditions/IsocohedronI.csv", T, err1);
    double T2;
    bool err2;
    matrix<double> pos = importcsv("/home/dino/Documents/tylercollab/Repo/Code/Basic/InitialConditions/IsocohedronP.csv", T2, err2);
    double k = 10.0;
    double rm = 1.25;

    B.k = k;
    B.rm = rm;
    B.par = pairs;
    B.posi = pos;

    string stringbase = "Num_mon=";
    stringstream ss;
    ss << monomers;
    string ss2 = ss.str();

    stringbase += ss2;

    A.run_with_real_surface(1000000,1000,B,stringbase);
    // A.run(1000000, 1000);

    return 0;
}