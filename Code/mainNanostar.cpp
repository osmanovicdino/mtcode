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
#include "MathHelper/Helper.h"
#include "MathHelper/Helper.cpp"
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

    matrix<double> K(3,3);

    matrix<double> K2 = K*K;

    cout << endl;
    // srand(time(NULL));


    Nanostar A(1,11.0);
    cout << "cons works" << '\n';

    cout << A.obj -> getdat() << endl;


    // A.Passa_set_nanostar(start, 30, 20, 4, 3, 5, "test.csv");
    //
    //
    // A.sortPairsTriplets(4, 3);
    //
    // for (int i = 0; i < A.bindpairs.size(); i++)
    // {
    //   cout << A.bindpairs[i] << '\n';
    // }
    // for (int i = 0; i < A.bendtriples.size(); i++)
    // {
    //   cout << A.bendtriples[i].a << A.bendtriples[i].b << A.bendtriples[i].c  << '\n';
    // }
    // # of timesteps,
    A.run(10000, 100);
    //
    // A.initStickerList(4, 3);
    // for (int i = 0; i<A.stickerList.size();i++)
    // {
    //   cout << A.stickerList[i] << '\n';
    // }
    //
    // vector<mdpair> testSet;
    // mdpair test1;
    // test1.a = 3;
    // test1.b = 6;
    // mdpair test2;
    // test2.a = 1;
    // test2.b = 2;
    // testSet.push_back(test1);
    // testSet.push_back(test2);
    //
    // vector<mdpair> output = A.inStickerList(testSet);
    // cout << output[0] << '\n';

    // A.run(10000,100);

    // string filename ="./";

    //int a = system("python3 /home/dino/Documents/Condensate/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

    //A.set_initial_state(filename);

    // for(int i = 0  ; i < A.bindpairs.size() ; i++) {
    // mdpair temp = A.bindpairs[i];
    // cout << temp.a << " " << temp.b << endl;
    // }

    // A.run(10000,100);


}
