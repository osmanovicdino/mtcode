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
#include "MDBase/Analysis/PairCorrelation.cpp"

// #include "NCGasR.h"
// #include "Microtubule.h"

//#include "MDGPU.cu"

using namespace std;

int main(int argc, char **argv)
{

    srand(time(NULL));


    //s_matrix<double> stat = importcsv("/home/dino/Code/DinoMD2/sim-19-11-29-16:55:47/x400_eps=0_kT=1.csv",T,err1);

    // string item_name;
    // ifstream nameFileout;
    // nameFileout.open("/home/dino/Code/DinoMD2/sim-19-12-20-17:05:04/files.txt");
    // string line;

    string preamable = "/media/dino/My";
    preamable+=" ";
    preamable+="Passport/Work_PC_MIT/Code/DinoFastMD/sim-20-04-16-13:08:06/";


 

    string temp = "a";

    vector1<string> a(5, temp);
    a[0] = "0";
    a[1] = "25";
    a[2] = "50";
    a[3] = "75";
    a[4] = "100";

    vector1<string> b(3, temp);
    b[0] = "8000";
    b[1] = "10000";
    b[2] = "12000";

    stringstream ss;

    int i = 9999;
    ss << setw(4) << setfill('0') << i;

    int i1 = 2;
    int i2 = 1;
    string xstring = preamable + "x" + ss.str() + "_v0=" + a[i1] + "_N=" + b[i2] + "_kT=1_l=140.csv";
    string oristring = preamable + "ori" + ss.str() + "_v0=" + a[i1] + "_N=" + b[i2] + "_kT=1_l=140.csv";

    double max = 5.0;
    double dx = 0.1;

    vector1<bool> pb(2,true);
    cube bc(140.,pb,2);

    double T;
    bool err1,err2;



    matrix<double> dat = importcsv(xstring,T,err1);

    matrix<double> ori = importcsv(oristring, T, err1);


    vector1<double> orig(ori.getNsafe());
    for(int k = 0 ; k < orig.getsize() ; k++) {
        orig[k] = ori(k,0);
    }

    matrix<int> res = calcg1(dat,orig,bc, max,dx,(2*pi/20.));

    outfunc(res,"res");

        // ofstream myfile;
        // myfile.open("/home/dino/Code/DinoMD2/sim-19-12-20-17:05:04/gs.csv");
        // while(std::getline(nameFileout, line))
        // {
        // s_matrix<double> stat = importcsv(line,T,err1);

        // vector1<int> g1 = calcg1(stat,bc,max,dx);
        // //cout << g1 << endl;
        // myfile <<= g1;
        // myfile << endl;
        // //vector1<double> pp = densitycorrelations(stat,bc,max,30);
        //    //std::cout << "line:" << line << std::endl;
        //    // TODO: assign item_name based on line (or if the entire line is
        //    // the item name, replace line with item_name in the code above)
        // }
        // myfile.close();

        //pausel();
        // s_matrix<double> stat = importcsv("/home/dino/Code/DinoMD2/sim-19-11-29-16:55:47/x400_eps=0_kT=1.csv",T,err1);

        // vector1<bool> pb(3,true);
        // cube bc(29.693,pb,3);
        // ofstream myfile;
        // myfile.open("/home/dino/Code/DinoMD2/sim-19-12-20-17:05:04/gs.csv");
        // while(std::getline(nameFileout, line))
        // {
        // s_matrix<double> stat = importcsv(line,T,err1);

        // vector1<int> g1 = calcg1(stat,bc,max,dx);
        // //cout << g1 << endl;
        // myfile <<= g1;
        // myfile << endl;

        return 0;
}