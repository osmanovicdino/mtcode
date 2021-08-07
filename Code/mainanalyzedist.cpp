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
#include <string.h>

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
#include "MDBase/Analysis/AnalysisFunctions.h"

using namespace std;

void morphNumericString(char *s, int n)
{
    char *p;
    int count;

    p = strchr(s, '.'); // Find decimal point, if any.
    if (p != NULL)
    {
        count = n; // Adjust for more or less decimals.
        while (count >= 0)
        { // Maximum decimals allowed.
            count--;
            if (*p == '\0') // If there's less than desired.
                break;
            p++; // Next character.
        }

        *p-- = '\0';      // Truncate string.
        while (*p == '0') // Remove trailing zeros.
            *p-- = '\0';

        if (*p == '.')
        { // If all decimals were zeros, remove ".".
            *p = '\0';
        }
    }

    cout << p << endl;
}

void add_double_to_ss(stringstream &ss, double d) {
    double intpart;
    if (modf(d, &intpart) == 0.0)
    {
        ss << d << ".";
    }
    else
    {
        std::string s;
        std::string t;
        std::stringstream out;
        out << d;
        s = out.str();

        t = s.substr(s.find(".") + 1);
        //cout << "number of decimal places: " << t.length();

        ss << std::setprecision(t.length()) << d;
    }
}

int main(int argc, char **argv)
{

srand(time(NULL));
signal(SIGSEGV, handler);
string mydir;
double l;
double bin;
//string mydir = "/u/scratch/d/dinoo/GrowthRun1/";

if (argc == 4)
{
    mydir = string(argv[1]);
    l = atof(argv[2]);
    bin = atof(argv[3]);
}
else{
    error("incorrect number of arguments for this executable");
}
    //mydir += base;//den=0.001_d=1._e=12._a=0.9272952180016123_arms=3";
    //cout << "all good" << endl;
    //"den=0.01_d=1._e=12._a=1.1592794807274085_arms=3"

    // string mydir2 = "den=0.001_d=1._e=12._a=1.1592794807274085_arms=3";

    
    matrix<int> growcurve = getgrowthcurve_distance_periodic_subset(mydir, l, bin, 5, 4000, 16000);

    string gc = "/growthwa";

    cout << mydir + gc << endl;
    cout << "output" << endl;

    outfunc(growcurve, mydir + gc);
}