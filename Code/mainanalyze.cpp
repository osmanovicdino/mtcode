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

int runtime;
double packing_fraction;
double size_of_part;
int num_arms;
double angs;
double energ;



if (argc == 7)
{
    packing_fraction = atof(argv[1]);
    size_of_part = atof(argv[2]);
    energ = atof(argv[3]);
    angs = atof(argv[4]);
    num_arms = atof(argv[5]);
}
else
{
    error("incorrect arg number");
}



// char buffer[50];
// int n;
// double a = 5, b = 3;
// n = sprintf(buffer, "%.20f", angs);
// printf("[%s] is a string %d chars long\n", buffer, n);
// morphNumericString(buffer,n);

// pausel();


int nt = 10000;
vector1<int> vec1(1);
vec1[0] = num_arms;

vector1<int> numb(1);

numb[0] = nt;

int tot = 0;
for (int i = 0; i < 1; i++)
{
    for (int j = i; j < 1; j++)
    {
        tot += vec1[i] * vec1[j];
    }
}

cout << tot << endl;

matrix<double> params2(tot * tot, 3);

for (int i = 0; i < tot * tot; i++)
{
    params2(i, 0) = energ;
    params2(i, 1) = 0.4 + size_of_part;
    params2(i, 2) = angs;
}

matrix<double> orient(num_arms, 3);

double nx1 = 0.0;
double ny1 = 1.;
double nz1 = 0.0;

double nx2 = -sqrt(12. / 16.);
double ny2 = -1. / 2.;
double nz2 = 0.0;

double nx3 = sqrt(12. / 16.);
double ny3 = -1. / 2.;
double nz3 = 0;

matrix<double> asd(3, 3);

asd(0, 0) = nx1;
asd(0, 1) = ny1;
asd(0, 2) = nz1;

asd(1, 0) = nx2;
asd(1, 1) = ny2;
asd(1, 2) = nz2;

asd(2, 0) = nx3;
asd(2, 1) = ny3;
asd(2, 2) = nz3;

double nx1_2 = sqrt(8. / 9.);
double ny1_2 = 0.;
double nz1_2 = -1. / 3.;

double nx2_2 = -sqrt(2. / 9.);
double ny2_2 = sqrt(2. / 3.);
double nz2_2 = -1. / 3.;

double nx3_2 = -sqrt(2. / 9.);
double ny3_2 = -sqrt(2. / 3.);
double nz3_2 = -1. / 3.;

double nx4_2 = 0;
double ny4_2 = 0;
double nz4_2 = 1.;

matrix<double> asd2(4, 3);

asd2(0, 0) = nx1_2;
asd2(0, 1) = ny1_2;
asd2(0, 2) = nz1_2;

asd2(1, 0) = nx2_2;
asd2(1, 1) = ny2_2;
asd2(1, 2) = nz2_2;

asd2(2, 0) = nx3_2;
asd2(2, 1) = ny3_2;
asd2(2, 2) = nz3_2;

asd2(3, 0) = nx4_2;
asd2(3, 1) = ny4_2;
asd2(3, 2) = nz4_2;

int iter2 = 0;
for (int i = 0; i < 1; i++)
{
    for (int j = 0; j < num_arms; j++)
    {
        if (num_arms == 3)
        {
            orient(iter2, 0) = asd(j, 0);
            orient(iter2, 1) = asd(j, 1);
            orient(iter2, 2) = asd(j, 2);
            iter2++;
        }
        else if (num_arms == 4)
        {
            orient(iter2, 0) = asd2(j, 0);
            orient(iter2, 1) = asd2(j, 1);
            orient(iter2, 2) = asd2(j, 2);
            iter2++;
        }
    }
}

GeneralPatch c4(vec1, numb, params2, orient);



string base = "den=";
stringstream dd;
dd << packing_fraction;
base += dd.str();

base += "_d=";
stringstream ss;
add_double_to_ss(ss,size_of_part);
base += ss.str();

stringstream ss1;
add_double_to_ss(ss1, energ);
base += "_e=";
base += ss1.str();

base += "_a=";
stringstream ss2;
add_double_to_ss(ss2, angs);
base += ss2.str();

base += "_arms=";
stringstream ii4;
ii4 << num_arms;
base += ii4.str();

string mydir = "/u/scratch/d/dinoo/GrowthRun1/";
//string mydir = "/home/dino/External/GrowthRun1/";
//mydir += base;//den=0.001_d=1._e=12._a=0.9272952180016123_arms=3";
//cout << "all good" << endl;


// string mydir2 = "den=0.001_d=1._e=12._a=1.1592794807274085_arms=3";
string mydir2 = string(argv[6]);

mydir += mydir2;
matrix<int> growcurve= getgrowthcurve(mydir, c4, 5);



string gc = "/growth";

cout << mydir + gc << endl;
cout << "output" << endl;


outfunc(growcurve,mydir+gc);

}