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
#include "Nanotube/Nanotube.h"

// #include "NCGasR.h"
// #include "Microtubule.h"

// #include "MDGPU.cu"

using namespace std;


int main(int argc, char **argv)
{
    unsigned long seed = mix(clock(), time(NULL), getpid());
    srand(seed);
    // int NM2;
    // double deltaG2,angle2;
    // A.run_with_real_surface(100000000, 10000, B, constantF, stringbase);
    // A.run(1000000, 1000);

    


    // s[0] = "2";
    // s[1] = "2110000";
    // s[2] = "2110000";
    // s[3] = "1111";
    // s[4] = "1111";
    // s[5] = "1111";
    // s[6] = "0";
    // s[7] = "11";
    string importstring;
    if (argc == 2)
    {
        stringstream ss;
        ss << argv[1];
        importstring = ss.str();
    }
    else
    {
        error("no");
    }

    // ofstream myfile;
    // myfile.open("g.csv");
    // for(int i = 0  ; i < s.size() ; i++) {
    //     myfile << s[i] << endl;
    // }
    // myfile.close();
    std::ifstream file(importstring); // your CSV file
    std::vector<std::string> s;
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Error: could not open file.\n";
        return 1;
    }

    while (std::getline(file, line))
    {
        s.push_back(line);
    }

    file.close();

    geneticcode g(s);

    // g.rate = 0.0;


    bool cc;
    NanotubeAssembly A(20., (g.no_types)*5000+1,cc);



    GeneralPatch c(CreateGeneralPatch(100., 1, 1.2, 0.6, g));


    
    A.setpots(c);
    A.setkT(1.0);
    A.setviscosity(1.0);

    A.run_box_equil(10000000, 1000, 100., g, "");
    // cout << a.no_types << endl;
    // cout << *(a.patch_num) << endl;
    // cout << *(a.patch_pos) << endl;
   
    // for(int j = 0  ; j < a.interactions[0].size() ; j++ )
    // cout << a.interactions[0][j] << ",";
    // cout << endl;

    // cout << a.rate << endl;
    // cout << *(a.proportions) << endl;



    // particle_adder vv;

    // box_vol vol2;
    // double ll = 20.;
    // vol2.ll = ll;
    // vol2.llx = ll;
    // vol2.lly = ll;
    // vol2.llz = 5.;
    // vv.set_volume(vol2);
    // cout << a.rate << endl;
    // vv.set_rate(a.rate);
    // vector1<int> weights((*(a.proportions)));
    // vv.set_irreducible_weights(weights);

    // vector<int> indices_to_add;
    // for(int i = 2 ; i < 20000 ; i++) {
    //     indices_to_add.push_back(i);
    // }
    // vector1<int> counts(a.no_types+1);
    // counts[0]=0;
    // counts[1]=4999;
    // counts[2]=9999;
    // counts[3]=14999;
    // counts[4]=19999;

    // vv.set_indices(indices_to_add);
    // vv.set_counts(counts);

    // double radius = 40;
    // // cout << "done 1" << endl;
    // bool cc;
    // NanotubeAssembly A(ll, 20000,cc);

    // // cout << "done 2" << endl;

    // vector<int> indices;
    // indices.push_back(0);
    // indices.push_back(1);

    
    // vector1<double> ve(3);
    // int fi;

    // cout << "add particles" << endl;
    // for(int i = 0  ; i < 10000 ; i++) {

    // bool does_return = false;
    // cout << vv.counts << endl;
    // vv.add_p_w(*(A.obj), indices, does_return, ve,fi);
    
    // if(does_return) {
    // cout << ve << endl;
    // cout << fi << endl;
    // cout << vv.counts << endl;

    
    // pausel();
    // }

    // }

    //we need to turn our genetic code into arguments
    
    return 0;
}