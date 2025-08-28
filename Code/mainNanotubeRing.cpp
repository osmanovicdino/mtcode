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
// #include <thrust/host_vector.h>
// #include <thrust/device_vector.h>
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
// #include "MDBase/potential.h"
#include "MDBase/MD.h"
#include "MDBase/Langevin.h"
// #include "MDBase/LangevinR.h"
#include "Condensate/Condensate.h"

// #include "NCGasR.h"
// #include "Microtubule.h"

// #include "MDGPU.cu"

using namespace std;

int main(int argc, char **argv)
{

    srand(time(NULL));



    double deltaG2, angle2;
    if (argc == 3)
    {
        deltaG2 = atof(argv[1]);
        angle2 = atof(argv[2]);
     }
    else
    {
        error("specify number of monomers");
    }
    
    int NM = 100;
    int NM2 = 0;
    // nm2 is the total amount of added crosslinks

    // signal(SIGSEGV, handler);

    Condensate A(10.,NM);


    // ShellProperties B;


    // stringstream sx1;
    // stringstream sx2;

    // sx1 << Ns;
    // sx2 << targetdensity;

    // RUN THIS TO REGENERATE SPHERE POINTS
    // string basic ="./Plotting/GenerateSpherePoints.wls";

    // string basic = "/home/dino/Documents/tylercollab/Repo/Code/Plotting/GenerateSpherePoints.wls";
    // string gap = " ";

    // string command = basic + gap + sx1.str() + gap + sx2.str();

    // system(command.c_str());

    // pausel();

    // int T;
    // bool err1;
    // matrix<int> pairs = importcsv("./IsocohedronI.csv", T, err1);
    // double T2;
    // bool err2;
    // matrix<double> pos = importcsv("./IsocohedronP.csv", T2, err2);
    // double k = 5.0;
    // double rm = 0.;

    // system("rm Iso*.csv");
    // B.k = k;
    // B.rm = rm;
    // B.par = pairs;
    // B.posi = pos;

    // B.DoAnMC(100.,false);

    // outfunc(B.posi,"res");
    // pausel();

    // vector1<double> mean = meanmat_end(pos,0);

    double approxradius = 20.0;

    double radius = 1.1 * 2 * approxradius;
    double monomers = NM;

    // monomers/(4/3piradius3)

    cout << "starting" << endl;
    // NanotubeAssembly A(radius, monomers);
    
    double deltaG = deltaG2;
    double angle = angle2;
    // BivalentPatch c2(deltaG, 1.4, angle);
    vector1<int> vec1(2);
    vec1[0] = 4;
    vec1[1] = 2;

    matrix<double> orient(6, 3);



    double nx1 = 1.0; //-0.5;
    double ny1 = 0.0; //-0.5 * sqrt(3.);
    double nz1 = 0.0; // 0.;

    double nx2 = -0.92388;  // 0.5;
    double ny2 = 0.382683;  // 0.5 * sqrt(3);
    double nz2 = 0.;  // 0.0;

    double nx3 = 0.19509;  // 0.5;
    double ny3 = 0.980785; // 0.5 * sqrt(3);
    double nz3 = 0.;  // 0.0;

    double nx4 = -0.19509;  // 0.5;
    double ny4 = -0.980785; // 0.5 * sqrt(3);
    double nz4 = 0.;  // 0.0;

    double nx5 = 1.0;
    double ny5 = 0.0;
    double nz5 = 0.0;

    double nx6 = -1.0;
    double ny6 = 0.0;
    double nz6 = 0.0;

    orient(0, 0) = nx1;
    orient(0, 1) = ny1;
    orient(0, 2) = nz1;

    orient(1, 0) = nx2;
    orient(1, 1) = ny2;
    orient(1, 2) = nz2;

    orient(2, 0) = nx3;
    orient(2, 1) = ny3;
    orient(2, 2) = nz3;

    orient(3, 0) = nx4;
    orient(3, 1) = ny4;
    orient(3, 2) = nz4;

    orient(4, 0) = nx5;
    orient(4, 1) = ny5;
    orient(4, 2) = nz5;

    orient(5, 0) = nx6;
    orient(5, 1) = ny6;
    orient(5, 2) = nz6;

    int tot = 4 * 4 + 4 * 2 + 2 * 2;
    matrix<double> params(tot, 3);
    double range = 1.2;
    vector1<int> torsion(tot);

    int iter = 0;
    for (int i = 0; i < 4; i++) // nanostar/nanostar interaction
    {
        for (int j = 0; j < 4; j++)
        {
            if (i <2 && j < 2 && i != j) // the sides cannot interact
            {
                params(iter, 0) = deltaG;
                params(iter, 1) = range;
                params(iter, 2) = angle;
                torsion[iter]=1;
                iter++;
            }
            // else if (i != j) // we want it to be directional
            // {
            //     params(iter, 0) = 0.0;
            //     params(iter, 1) = range;
            //     params(iter, 2) = angle;
            //     iter++;
            // }
            else
            {
                params(iter, 0) = 0;
                params(iter, 1) = range;
                params(iter, 2) = angle;
                iter++;
            }
        }
    }

    for (int i = 0; i < 4; i++) // crosslinks
    {
        for (int j = 0; j < 2; j++)
        {
            if (i > 2)
            {
                params(iter, 0) = deltaG;
                params(iter, 1) = range;
                params(iter, 2) = angle;
            }
            else
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = range;
                params(iter, 2) = angle;
            }
            iter++;
        }
    }

    for (int i = 0; i < 2; i++) // caps
    {
        for (int j = 0; j < 2; j++)
        {
            if (i == 0 && j == 0) // only binds to one end
            {
                params(iter, 0) = 0.0;
                params(iter, 1) = range;
                params(iter, 2) = angle;
            }
            else
            {
                params(iter, 0) = 0.0; // bind to one end, blocking further growth
                params(iter, 1) = range;
                params(iter, 2) = angle;
            }
            iter++;
        }
    }


    vector1<int> numb(2);
    numb[0] = NM ;
    numb[1] = NM + NM2 ;

    //    TetrahedralWithBivalent c2(params, Ns+NM2 , Ns + NM,orient,orient2); //set the difference to be  greater
    
    
    GeneralPatchTorsion c2(vec1, numb, params, orient, torsion);
    // c2.v = orient;
    // c2.v2= orient2;

    A.setpots(c2);
    // c2.v = orient;
    // c2.v2= orient2;

    // A.setpots(c2);
    A.setkT(1.0);

    cout << "done" << endl;

    // outfunc(B.posi,"res");
    // cout << "MC done" << endl;

    // pausel();

    string stringbase = "Num_mon=";
    stringstream ss;
    ss << monomers;
    string str = ss.str();

    stringbase += str;

    stringstream ss2;
    ss2 << deltaG;

    stringbase += string("dG=") + ss2.str();

    stringstream ss3;
    ss3 << angle;

    stringbase += string("ang=") + ss3.str();

    // matrix<double> constantF(Ns+NM,3);
    // constantF(0,2) = -100.;
    // constantF(4095,2) = 100.;
    // A.conf.setv(0.0); // no confinement
    //A.run(10000000, 10000, stringbase);


    A.run(400000,100,stringbase);
    // A.run_add_particles(10000000, 10000, 0.001, stringbase);
    //  A.run_with_real_surface(100000000, 10000, B, constantF, stringbase);
    //  A.run(1000000, 1000);

    return 0;
}