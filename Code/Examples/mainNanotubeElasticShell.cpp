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

    if (argc == 2)
    {
        NM = atof(argv[1]);
    }
    else
    {
        error("specify number of monomers");
    }

    // signal(SIGSEGV, handler);

    ShellProperties B;
    int Ns = 2048;
    double targetdensity = 2.0;

    stringstream sx1;
    stringstream sx2;

    sx1 << Ns;
    sx2 << targetdensity;

    // string basic ="./Plotting/GenerateSpherePoints.wls";

    string basic = "/home/dino/Documents/tylercollab/Repo/Code/Plotting/GenerateSpherePoints.wls";
    string gap = " ";

    string command = basic + gap + sx1.str() + gap + sx2.str();

    system(command.c_str());

    int T;
    bool err1;
    matrix<int> pairs = importcsv("./IsocohedronI.csv", T, err1);
    double T2;
    bool err2;
    matrix<double> pos = importcsv("./IsocohedronP.csv", T2, err2);
    double k = 5.0;
    double rm = 1.25;

    system("rm Iso*.csv");
    B.k = k;
    B.rm = rm;
    B.par = pairs;
    B.posi = pos;

    // B.DoAnMC(100.,false);

    // outfunc(B.posi,"res");
    // pausel();

    // vector1<double> mean = meanmat_end(pos,0);

    double approxradius = sqrt(SQR(pos(0, 0)) + SQR(pos(0, 1)) + SQR(pos(0, 2)));

    double radius = (1. / 0.9) * 2 * approxradius;
    double monomers = NM;

    // monomers/(4/3piradius3)

    cout << "starting" << endl;
    NanotubeAssembly A(radius, monomers);

    double deltaG = 30.0;
    double angle = 0.9;
    // BivalentPatch c2(deltaG, 1.4, angle);

    matrix<double> orient(5, 3);

    double nx4 = 1.0;
    double ny4 = 0.0;
    double nz4 = 0.0;

    double nx5 = -1.0;
    double ny5 = 0.0;
    double nz5 = 0.0;

    double nx6 = -0.5;
    double ny6 = 0.5 * sqrt(3.);
    double nz6 = 0.;

    double nx7 = -0.5;
    double ny7 = -0.5 * sqrt(3.);
    double nz7 = 0.;

    double nx8 = 0.5;
    double ny8 = 0.5 * sqrt(3);
    double nz8 = 0.0;

    orient(0, 0) = nx6;
    orient(0, 1) = ny6;
    orient(0, 2) = nz6;

    orient(1, 0) = nx7;
    orient(1, 1) = ny7;
    orient(1, 2) = nz7;

    orient(2, 0) = nx8;
    orient(2, 1) = ny8;
    orient(2, 2) = nz8;

    orient(3, 0) = nx4;
    orient(3, 1) = ny4;
    orient(3, 2) = nz4;

    orient(4, 0) = nx5;
    orient(4, 1) = ny5;
    orient(4, 2) = nz5;

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

    TetrahedralWithBivalent c2(params, Ns + NM / 4, Ns + NM);

    A.setpots(c2);
    A.setkT(1.0);

    cout << "done" << endl;

    // outfunc(B.posi,"res");
    // cout << "MC done" << endl;

    // pausel();

    string stringbase = "Num_mon=";
    stringstream ss;
    ss << monomers;
    string ss2 = ss.str();

    stringbase += ss2;

    A.run_with_real_surface_add_particles(10000000, 1000, B, 0.000, stringbase);
    // A.run(1000000, 1000);

    return 0;
}