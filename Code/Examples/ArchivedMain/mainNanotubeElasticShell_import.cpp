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
    
    
    string paramfile;
    string olddatfile;
    string oldorifile;
    string oldindfile;
    string shellpairsfile;
    if (argc == 5)
    {
        stringstream ss1,ss2,ss3,ss4,ss5;
        ss1 << argv[1];
        ss2 << argv[2];
        ss3 << argv[3];
        ss4 << argv[4];
        //ss5 << argv[5];

        paramfile = ss1.str();
        olddatfile = ss2.str();
        oldorifile = ss3.str();
        oldindfile = ss4.str();
        //shellpairsfile = ss5.str();
        // NM2 = atof(argv[1]);
        // deltaG2 = atof(argv[2]);
        // angle2 = atof(argv[3]);
    }
    else
    {
        error("specify number of monomers");
    }
    double T;
    bool err1;
    matrix<double> sim_params = importcsv(paramfile, T, err1);

    bool err2;
    matrix<double> olddat = importcsv(olddatfile, T, err2);

    bool erro;
    matrix<double> oldori = importcsv(oldorifile, T, erro);

    int T2;
    bool err3;
    matrix<double> oldind_temp = importcsv(oldindfile, T2, err3);

    int Td;
    bool err4;
    matrix<int> pairs = importcsv("./IsocohedronI.csv", Td, err4);

    int Ti;
    bool erri;
    matrix<int> quads = importcsv("./IsocohedronI2.csv", Ti, erri);

    double k = sim_params(1, 0);

    double rm = sim_params(1, 1);

    double kappa = sim_params(1, 3);

    double Td4;
    bool errd4;
    matrix<double> bindingdis = importcsv("./IsocohedronD.csv", Td4, errd4);

    ShellProperties B(pairs, quads, olddat, bindingdis, k, rm, kappa);

    if(err1 || err2 || err3 || err4 || erri || errd4 ) {
        cout << paramfile << " " << err1 << endl;
        cout << olddatfile << " " << err2 << endl;
        cout << oldorifile << " " << err3 << endl;
        cout << oldindfile << " " << err4 << endl;
        cout << shellpairsfile << " " << erro << endl;

        error("files not imported correctly");

    }

    int NM = sim_params(0, 0);
    //nm2 is the total amount of added crosslinks

    // signal(SIGSEGV, handler);
    vector1<int> oldind(oldind_temp.getnrows());
    for (int i = 0; i < oldind_temp.getnrows(); i++)
    {
        oldind[i] = oldind_temp(i, 0);
    }



    // ShellProperties B;
    int Ns = sim_params(0, 1);
    //double targetdensity = 2.0;
    if(olddat.getnrows() < Ns) error("size of shell larger than number of monomers being imported");
    if(oldind_temp.getnrows() < Ns) error("size of index file is mismatched");
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



    // double T2;
    // bool err2;
    // matrix<double> pos = importcsv("./IsocohedronP.csv", T2, err2);


    matrix<double> pos(Ns,3);
    for(int i = 0  ; i < Ns ; i++) {
        for(int j = 0  ; j < 3 ; j++) {
            pos(i,j)=olddat(i,j);
        }
    }

    // system("rm Iso*.csv");
    B.k = k;
    B.rm = rm;
    B.par = pairs;
    B.posi = pos;

    // B.DoAnMC(100.,false);

    // outfunc(B.posi,"res");
    // pausel();

    // vector1<double> mean = meanmat_end(pos,0);
    
//    double approxradius = sqrt(SQR(pos(0, 0)) + SQR(pos(0, 1)) + SQR(pos(0, 2)));

    double radius = sim_params(1, 2);
    double monomers = NM;

    // monomers/(4/3piradius3)

    NanotubeAssembly A(radius, monomers);
    NM2 = sim_params(2,0);
    int NM3=sim_params(2,1);
    double deltaG = sim_params(2, 2);
    double angle = sim_params(2, 3);
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
            if (i == 2 || j == 2)
            {
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
            if (i == 2)
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
            params(iter, 0) = 0.0; // deltaG;
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
    vector<string> posfiles;
    cout << "ok1" << endl;
    return_csv_in_current_dir("pos",posfiles);
    cout << "ok2" << endl;
    double prod = sim_params(3, 0);
    WeiM c1;
    c1.M = sim_params(3, 1);
    c1.weight = sim_params(3, 2);
    A.run_with_real_surface_add_particles_continue(10000000, 10000, posfiles.size(), B, prod, c1, olddat,oldori,oldind,stringbase);
    // A.run_with_real_surface(100000000, 10000, B, constantF, stringbase);
        // A.run(1000000, 1000);

        return 0;
}