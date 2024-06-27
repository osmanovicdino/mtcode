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
    // int NM2;
    // double deltaG2,angle2;
    string importstring;
    
    if (argc == 2)
    {
        stringstream ss;
        ss << argv[1];
        importstring = ss.str();


    }
    else
    {
        error("specify input file");
    }

    double T3;
    bool err3;
    matrix<double> sim_params = importcsv(importstring, T3, err3);

    if(err3) error("param file not imported correctly");

    int NM = sim_params(0,0);
    //nm2 is the total amount of added crosslinks

    // signal(SIGSEGV, handler);

    
    int Ns = 4322; //technically this should be overwrittable, however, because we generated the circle with mathematica we leave it fixed
    double targetdensity = 2.0;

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
    matrix<int> quads = importcsv("./IsocohedronI2.csv", Ti, erri);
    double T2;
    bool err2;
    matrix<double> pos = importcsv("./IsocohedronP2.csv", T2, err2);
    double k = sim_params(1, 0);
    
    double rm = sim_params(1, 1);

    double kappa = sim_params(1, 3);

    double T4;
    bool err4;
    matrix<double> bindingdis = importcsv("./IsocohedronD.csv", T4, err4);

    ShellProperties B(pairs,quads,pos,bindingdis,k,rm,kappa);

    // system("rm Iso*.csv");
    // B.k = k;
    // B.rm = rm;
    // B.par = pairs;
    // B.posi = pos;

    // B.DoAnMC(100.,false);

    // outfunc(B.posi,"res");
    // pausel();

    // vector1<double> mean = meanmat_end(pos,0);
    
    // double approxradius = sqrt(SQR(pos(0, 0)) + SQR(pos(0, 1)) + SQR(pos(0, 2)));


    double radius = sim_params(1,2);
    double monomers = NM;

    // monomers/(4/3piradius3)

    cout << "starting" << endl;
    NanotubeAssembly A(radius, monomers);
    int NM2 =  sim_params(2,0);
    int NM3 = sim_params(2, 1);
    double deltaG = sim_params(2,2);
    double angle = sim_params(2,3);
    // BivalentPatch c2(deltaG, 1.4, angle);
    vector1<int> vec1(3);
    vec1[0] = 3;
    vec1[1] = 2;
    vec1[2] = 2;

    matrix<double> orient(7, 3);

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

    orient(3, 0) = nx4;
    orient(3, 1) = ny4;
    orient(3, 2) = nz4;

    orient(4, 0) = nx5;
    orient(4, 1) = ny5;
    orient(4, 2) = nz5;

    orient(5, 0) = nx4;
    orient(5, 1) = ny4;
    orient(5, 2) = nz4;

    orient(6, 0) = nx5;
    orient(6, 1) = ny5;
    orient(6, 2) = nz5;


    int tot = 3 * 3 + 3 * 2 + 3 * 2 + 2 * 2 + 2 * 2 + 2 * 2;
    matrix<double> params(tot, 3);
    double range = 1.2;
  
    int iter = 0;
    for (int i = 0; i < 3; i++) // nanostar/nanostar interaction
    {
        for (int j = 0; j < 3; j++)
        {
            if(i==2||j==2) {
                params(iter, 0) = 0.0;
                params(iter, 1) = range;
                params(iter, 2) = angle;
                iter++;
            }
            else 
            {
                params(iter, 0) = deltaG;
                params(iter, 1) = range;
                params(iter, 2) = angle;
                iter++;
            }
        }
    }

    for (int i = 0; i < 3; i++) // nanostar/anti-invader interaction
    {
        for (int j = 0; j < 2; j++)
        {
            if(i==2) {
            params(iter, 0) = deltaG;
            params(iter, 1) = range;
            params(iter, 2) = angle;
            }
            else{
            params(iter, 0) = 0.0;
            params(iter, 1) = range;
            params(iter, 2) = angle;
            }
            iter++;
        }
    }

    for (int i = 0; i < 3; i++) // nanonstar/invader interaction
    {
        for (int j = 0; j < 2; j++)
        {
            if (i == 2)
            {
            params(iter, 0) = 0.;
            params(iter, 1) = range;
            params(iter, 2) = angle;
            }
            else
            {
            params(iter, 0) = deltaG;
            params(iter, 1) = range;
            params(iter, 2) = angle;
            }
            iter++;
        }
    }
    for (int i = 0; i < 2; i++) // anti-invader/anti-invader interaction
    {
        for (int j = 0; j < 2; j++)
        {
            params(iter, 0) = 0.0;//deltaG;
            params(iter, 1) = range;
            params(iter, 2) = angle;

            iter++;
        }
    }

    for (int i = 0; i < 2; i++) // anti-invader/invader interaction
    {
        for (int j = 0; j < 2; j++)
        {
            params(iter, 0) = 0.0;
            params(iter, 1) = range;
            params(iter, 2) = angle;

            iter++;
        }
    }

    for (int i = 0; i < 2; i++) // invader/invader interaction
    {
        for (int j = 0; j < 2; j++)
        {
            params(iter, 0) = 0.0;//deltaG
            params(iter, 1) = range;
            params(iter, 2) = angle;

            iter++;
        }
    }
    vector1<int> numb(3);
    numb[0]=Ns+NM2;
    numb[1]=Ns+NM2+NM3;
    numb[2]=Ns+NM;

 
 
//    TetrahedralWithBivalent c2(params, Ns+NM2 , Ns + NM,orient,orient2); //set the difference to be  greater
    GeneralPatch c2(vec1, numb, params, orient);
    // c2.v = orient;
    // c2.v2= orient2;


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
    // matrix<double> constantF(Ns+NM,3);
    // constantF(0,2) = -100.;
    // constantF(4095,2) = 100.;
    double prod = sim_params(3,0);
    WeiM c1;
    c1.M=sim_params(3,1);
    c1.weight=sim_params(3,2);

    A.run_with_real_surface_add_particles(10000000, 10000, B, prod, c1, stringbase);
    // A.run_with_real_surface(100000000, 10000, B, constantF, stringbase);
        // A.run(1000000, 1000);

        return 0;
    }