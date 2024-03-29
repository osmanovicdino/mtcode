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

    srand(time(NULL));

    int arm_length;
    if(argc == 2) {
        arm_length =  atof(argv[1]);
    }
    Nanostar A(1000,200.0,arm_length);


    //A.DoAnMC();

    vector<string> posfiles;
    return_csv_in_current_dir("pos", posfiles);
    int s2 = posfiles.size();
    double T;
    int TT;
    bool vv1, vv2, vv3;
    matrix<double> postemp = importcsv(posfiles[posfiles.size() - 1], T, vv1);

    cout << "everything seems to work" << endl;
    A.obj->setdat(postemp);
    //A.Passa_set_nanostar();

    A.run(1000000,1000,s2);

    

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