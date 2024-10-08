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
    int NM2;
    double deltaG2,angle2;
    double F1;
    double k2,rm2,kapp2;
    if (argc == 5)
    {
        F1 = atof(argv[1]);
        k2 = atof(argv[2]);
        rm2 = atof(argv[3]);
        kapp2 = atof(argv[4]);
    }
    else
    {
        error("specify number of monomers");
    }
    int NM = 10;
    //nm2 is the total amount of added crosslinks

    // signal(SIGSEGV, handler);
    NM2 = 2000;
    deltaG2 = 100.;
    angle2 = 0.6;


    //double targetdensity = 2.0;

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


    int T;
    bool err1;
    matrix<int> pairs = importcsv("./IsocohedronI.csv", T, err1);

    int Ti;
    bool erri;
    matrix<int> quads= importcsv("./IsocohedronI2.csv", Ti, erri);
    double T2;
    bool err2;
    matrix<double> pos = importcsv("./IsocohedronP2.csv", T2, err2);
    double k = k2;
    double rm = rm2;
    double kappa = kapp2;
    double T3;
    bool err3;
    matrix<double> bindingdis = importcsv("./IsocohedronD.csv", T3, err3);

    // system("rm Iso*.csv");

    ShellProperties B(pairs,quads,pos,bindingdis,k,rm,kappa);
    int Ns = pos.getnrows();


    // B.k = k;
    // B.rm = rm;
    // B.par = pairs;
    // B.posi = pos;
    // B.quad = quads;
    // B.kappa = kappa;

    // B.DoAnMC(100.,false);

    // outfunc(B.posi,"res");
    // pausel();

    // vector1<double> mean = meanmat_end(pos,0);
    
    double approxradius = 18.0;


    double radius = 1.1 * 2 * approxradius;
    double monomers = NM;

    // monomers/(4/3piradius3)

    cout << "starting" << endl;
    NanotubeAssembly A(radius, monomers);

    double deltaG = deltaG2;
    double angle = angle2;
    // BivalentPatch c2(deltaG, 1.4, angle);

    matrix<double> orient(4, 3);
    matrix<double> orient2(2, 3);

    double nx4 = 1.0;
    double ny4 = 0.0;
    double nz4 = 0.0;

    double nx5 = -1.0;
    double ny5 = 0.0;
    double nz5 = 0.0;

    double nx6 = 1.0;//-0.5;
    double ny6 = 0.0;//0.5 * sqrt(3.);
    double nz6 = 0.0;

    double nx7 = 0.0;//-0.5;
    double ny7 = 1.0;//-0.5 * sqrt(3.);
    double nz7 = 0.0;//0.;

    double nx8 = -1.;//0.5;
    double ny8 = 0.;//0.5 * sqrt(3);
    double nz8 = 0.;//0.0;

    double nx9 = 0.;//0.5;
    double ny9 = -1.;//0.5 * sqrt(3);
    double nz9 = 0.;//0.0;

    orient(0, 0) = nx6;
    orient(0, 1) = ny6;
    orient(0, 2) = nz6;

    orient(1, 0) = nx8;
    orient(1, 1) = ny8;
    orient(1, 2) = nz8;

    orient(2, 0) = nx7;
    orient(2, 1) = ny7;
    orient(2, 2) = nz7;

    orient(3, 0) = nx9;
    orient(3, 1) = ny9;
    orient(3, 2) = nz9;

    orient2(0, 0) = nx4;
    orient2(0, 1) = ny4;
    orient2(0, 2) = nz4;

    orient2(1, 0) = nx5;
    orient2(1, 1) = ny5;
    orient2(1, 2) = nz5;

    int tot = 4 * 4 + 4 * 2 + 2 * 2;
    int iter = 0;
    matrix<double> params(tot, 3);
    for (int i = 0; i < 4; i++)
    {
        for(int j = 0  ; j < 4 ; j++) {
        if(i<=1 && j <=1) {
        params(iter, 0) = 00.0;
        }
        else if (i >= 2 && j >= 2)
        {
        params(iter, 0) = 00.0;
        }
        else{
        params(iter, 0) = 00.0;
        }
        params(iter, 1) = 1.2;
        params(iter, 2) = angle;
        iter++;
        }
    }
    for (int i = 4 * 4; i < 4 * 4+4*2; i++)
    {
        if(i!=19 || i!=23) {
        params(i, 0) = deltaG; //interaction of 4s with 2
        params(i, 1) = 1.2;
        params(i, 2) = angle;
        }
        else{
        params(i, 0) = 0.0; // interaction of 4s with 2
        params(i, 1) = 1.2;
        params(i, 2) = angle;
        }
    }
    for (int i = 4 * 4 + 4 * 2; i < tot; i++)
    {
        params(i, 0) = deltaG; //interaction of backbone
        params(i, 1) = 1.2;
        params(i, 2) = angle;
    }

    TetrahedralWithBivalent c2(params, NM2 , Ns + NM,orient,orient2); //set the difference to be  greater

    // c2.v = orient;
    // c2.v2= orient2;


    A.setpots(c2);
    A.setkT(1.0);

    cout << "done" << endl;

    // outfunc(B.posi,"res");
    // cout << "MC done" << endl;

    // pausel();


    // stringstream ss2;
    // ss2 << deltaG;

    // stringbase += string("dG=")+ss2.str();

    // stringstream ss3;
    // ss3 << angle;
    string stringbase = "F=";
    stringstream ss;
    ss << 0.;
    string str = ss.str();

    stringstream ss2, ss3, ss4;
    ss2 << k;
    ss3 << rm;
    ss4 << kappa;

    stringbase += str;
    stringbase += string("_k=") + ss2.str();
    stringbase += string("_rm=") + ss3.str();
    stringbase += string("_kappa=") + ss4.str();

    matrix<double> constantF(Ns + NM, 3);

    A.run_with_real_surface(500000, 10000, B, constantF, stringbase, false);
    // stringbase += string("ang=")+ss3.str();
    for(double g=0 ; g <51. ; g+=10. ) {
    
    constantF(2,2) +=  g;
    constantF(5,2) -= g;

    string stringbase = "F=";
    stringstream ss;
    ss << F1+g;
    string str = ss.str();

    stringstream ss2, ss3, ss4;
    ss2 << k;
    ss3 << rm;
    ss4 << kappa;

    stringbase += str;
    stringbase += string("_k=") + ss2.str();
    stringbase += string("_rm=") + ss3.str();
    stringbase += string("_kappa=") + ss4.str();

    //A.conf.setv(0.0);//no confinement
    // A.run(20000000, 10000, stringbase);
        // A.run_add_particles(10000000, 10000, 0.001, stringbase);
         A.run_with_real_surface(500000, 10000, B, constantF, stringbase,true);
        //  A.run(1000000, 1000);
    }
        return 0;
}