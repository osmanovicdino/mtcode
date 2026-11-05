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
#include <unordered_set>
#include <algorithm>
#include <parallel/algorithm>
#include <array>
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

    vector<string> s;

    int str1 = 1 + rand() % 3;

    s.push_back(to_string(str1));
    vector1<int> opt1(6);
    opt1[0] = 1;
    vector1<int> opt2_1(6);
    opt2_1[0] = 1;
    opt2_1[1] = 1;
    vector1<int> opt2_2(6);
    opt2_2[0] = 1;
    opt2_2[2] = 1;

    vector1<int> opt3_1(6);
    vector1<int> opt3_2(6);
    opt3_1[0] = 1;
    opt3_1[1] = 1;
    opt3_1[2] = 1;

    opt3_2[0] = 1;
    opt3_2[2] = 1;
    opt3_2[4] = 1;

    vector1<int> opt4_1(6);
    vector1<int> opt4_2(6);
    opt4_1[0] = 1;
    opt4_1[1] = 1;
    opt4_1[2] = 1;
    opt4_1[3] = 1;

    opt4_2[0] = 1;
    opt4_2[1] = 1;
    opt4_2[2] = 1;
    opt4_2[4] = 1;

    vector1<int> opt5(6, 1);
    opt5[5] = 0;

    vector1<int> opt6(6, 1);

    // cout << opt1 << endl;

    std::vector<int> weights = {1, 4, 3, 2, 1, 1};

    // random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // discrete distribution with given weights
    std::discrete_distribution<> dist(weights.begin(), weights.end());

    matrix<int> olg(str1, 7);
    for (int i = 0; i < str1; i++)
    {

        int n = dist(gen) + 1;
        olg(i, 0) = n; // between 1 and 6

        // std::vector<int> v = {1, 2, 3, 4, 5, 6};

        // Pick random n between 2 and 6
        // std::random_device rd;
        // std::mt19937 gen(rd());
        // // std::uniform_int_distribution<> distN(2, v.size());
        // // int n = distN(gen);

        // // Shuffle and take first n
        // std::shuffle(v.begin(), v.end(), gen);

        // std::vector<int> chosen(v.begin(), v.begin() + n);
        // olg(i,1) = 1;
        // olg(i,2) = 1;
        // // for(int j = 0 ; j < chosen.size() ; j++ ) {
        // //     olg(i,chosen[j])=1;
        // // }
        // if(n==3) {
        //     olg(i,3) = 1; //identical to rotation
        // }
        // if(n==4) { //they can either all be on the same plane or out of plane
        //     int ra = rand() %2;
        //     if(ra == 0) {
        //         olg(i,3) = 1;
        //         olg(i,4) = 1;
        //     }
        //     else {
        //         olg(i,3) = 1;
        //         olg(i,5) = 1;
        //     }
        // }

        cout << n << endl;

        switch (n)
        {
        case 1:
        {
            for (int j = 0; j < 6; j++)
                olg(i, j + 1) = opt1[j];

            break;
        }
        case 2:
        {
            int randr = rand() % 2;
            if (randr == 0)
            {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt2_1[j];
            }
            else
            {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt2_2[j];
            }

            break;
        }
        case 3:
        {
            int randr = rand() % 2;
            if (randr == 0)
            {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt3_1[j];
            }
            else
            {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt3_2[j];
            }

            break;
        }
        case 4:
        {
            int randr = rand() % 2;
            if (randr == 0)
            {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt4_1[j];
            }
            else
            {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt4_2[j];
            }

            break;
        }
        case 5:
        {
            for (int j = 0; j < 6; j++)
                olg(i, j + 1) = opt5[j];

            break;
        }
        case 6:
        {
            for (int j = 0; j < 6; j++)
                olg(i, j + 1) = opt6[j];

            break;
        }
        }

        string result;

        result += std::to_string(n);
        for (int j = 1; j < 7; j++)
        {
            result += std::to_string(olg(i, j));
        }
        s.push_back(result);
    }

    for (int i = 0; i < str1; i++)
    {
        for (int j = i; j < str1; j++)
        {
            int totn = (olg(i, 0) * olg(j, 0));
            vector1<int> ints(totn, 1);
            for (int k = 0; k < totn; k++)
            {
                int ran = rand() % 6;
                if (ran == 5)
                    ints[k] = 0;
            }
            string result;

            for (int k = 0; k < totn; k++)
            {
                result += to_string(ints[k]);
            }
            s.push_back(result);
        }
    }

    s.push_back("-7");

    vector1<int> weig(str1);

    for (int i = 0; i < str1; i++)
    {
        weig[i] = 1 + rand() % 3;
    }

    string ress;

    for (int i = 0; i < str1; i++)
    {
        ress += to_string(weig[i]);
    }

    s.push_back(ress);

    // s[0] = "2";
    // s[1] = "2110000";
    // s[2] = "2110000";
    // s[3] = "1111";
    // s[4] = "1111";
    // s[5] = "1111";
    // s[6] = "0";
    // s[7] = "11";

    ofstream myfile;
    myfile.open("g.csv");
    for (int i = 0; i < s.size(); i++)
    {
        myfile << s[i] << endl;
    }
    myfile.close();

    geneticcode g(s);

    double ll =  20.;
    int Nfan = 229;

    std::array<std::array<double, 3>, 229> fan_data = {{{-5 + ll / 2, ll / 2, 11}, {-5 + ll / 2, ll / 2, 12}, {-5 + ll / 2, ll / 2, 13}, {-5 + ll / 2, ll / 2, 14}, {-5 + ll / 2, ll / 2, 15}, {-5 + ll / 2, ll / 2, 16}, {-5 + ll / 2, ll / 2, 17}, {-5 + ll / 2, ll / 2, 18}, {-5 + ll / 2, ll / 2, 19}, {-4 + ll / 2, ll / 2, 11}, {-4 + ll / 2, ll / 2, 12}, {-4 + ll / 2, ll / 2, 13}, {-4 + ll / 2, ll / 2, 14}, {-4 + ll / 2, ll / 2, 15}, {-4 + ll / 2, ll / 2, 16}, {-4 + ll / 2, ll / 2, 17}, {-4 + ll / 2, ll / 2, 18}, {-4 + ll / 2, ll / 2, 19}, {-3 + ll / 2, ll / 2, 11}, {-3 + ll / 2, ll / 2, 12}, {-3 + ll / 2, ll / 2, 13}, {-3 + ll / 2, ll / 2, 14}, {-3 + ll / 2, ll / 2, 15}, {-3 + ll / 2, ll / 2, 16}, {-3 + ll / 2, ll / 2, 17}, {-3 + ll / 2, ll / 2, 18}, {-3 + ll / 2, ll / 2, 19}, {-2 + ll / 2, ll / 2, 11}, {-2 + ll / 2, ll / 2, 12}, {-2 + ll / 2, ll / 2, 13}, {-2 + ll / 2, ll / 2, 14}, {-2 + ll / 2, ll / 2, 15}, {-2 + ll / 2, ll / 2, 16}, {-2 + ll / 2, ll / 2, 17}, {-2 + ll / 2, ll / 2, 18}, {-2 + ll / 2, ll / 2, 19}, {-1 + ll / 2, ll / 2, 11}, {-1 + ll / 2, ll / 2, 12}, {-1 + ll / 2, ll / 2, 13}, {-1 + ll / 2, ll / 2, 14}, {-1 + ll / 2, ll / 2, 15}, {-1 + ll / 2, ll / 2, 16}, {-1 + ll / 2, ll / 2, 17}, {-1 + ll / 2, ll / 2, 18}, {-1 + ll / 2, ll / 2, 19}, {ll / 2, ll / 2, 11}, {ll / 2, ll / 2, 12}, {ll / 2, ll / 2, 13}, {ll / 2, ll / 2, 14}, {ll / 2, ll / 2, 15}, {ll / 2, ll / 2, 16}, {ll / 2, ll / 2, 17}, {ll / 2, ll / 2, 18}, {ll / 2, ll / 2, 19}, {1 + ll / 2, ll / 2, 11}, {1 + ll / 2, ll / 2, 12}, {1 + ll / 2, ll / 2, 13}, {1 + ll / 2, ll / 2, 14}, {1 + ll / 2, ll / 2, 15}, {1 + ll / 2, ll / 2, 16}, {1 + ll / 2, ll / 2, 17}, {1 + ll / 2, ll / 2, 18}, {1 + ll / 2, ll / 2, 19}, {2 + ll / 2, ll / 2, 11}, {2 + ll / 2, ll / 2, 12}, {2 + ll / 2, ll / 2, 13}, {2 + ll / 2, ll / 2, 14}, {2 + ll / 2, ll / 2, 15}, {2 + ll / 2, ll / 2, 16}, {2 + ll / 2, ll / 2, 17}, {2 + ll / 2, ll / 2, 18}, {2 + ll / 2, ll / 2, 19}, {3 + ll / 2, ll / 2, 11}, {3 + ll / 2, ll / 2, 12}, {3 + ll / 2, ll / 2, 13}, {3 + ll / 2, ll / 2, 14}, {3 + ll / 2, ll / 2, 15}, {3 + ll / 2, ll / 2, 16}, {3 + ll / 2, ll / 2, 17}, {3 + ll / 2, ll / 2, 18}, {3 + ll / 2, ll / 2, 19}, {4 + ll / 2, ll / 2, 11}, {4 + ll / 2, ll / 2, 12}, {4 + ll / 2, ll / 2, 13}, {4 + ll / 2, ll / 2, 14}, {4 + ll / 2, ll / 2, 15}, {4 + ll / 2, ll / 2, 16}, {4 + ll / 2, ll / 2, 17}, {4 + ll / 2, ll / 2, 18}, {4 + ll / 2, ll / 2, 19}, {5 + ll / 2, ll / 2, 11}, {5 + ll / 2, ll / 2, 12}, {5 + ll / 2, ll / 2, 13}, {5 + ll / 2, ll / 2, 14}, {5 + ll / 2, ll / 2, 15}, {5 + ll / 2, ll / 2, 16}, {5 + ll / 2, ll / 2, 17}, {5 + ll / 2, ll / 2, 18}, {5 + ll / 2, ll / 2, 19}, {ll / 2, -5 + ll / 2, 11}, {ll / 2, -5 + ll / 2, 12}, {ll / 2, -5 + ll / 2, 13}, {ll / 2, -5 + ll / 2, 14}, {ll / 2, -5 + ll / 2, 15}, {ll / 2, -5 + ll / 2, 16}, {ll / 2, -5 + ll / 2, 17}, {ll / 2, -5 + ll / 2, 18}, {ll / 2, -5 + ll / 2, 19}, {ll / 2, -4 + ll / 2, 11}, {ll / 2, -4 + ll / 2, 12}, {ll / 2, -4 + ll / 2, 13}, {ll / 2, -4 + ll / 2, 14}, {ll / 2, -4 + ll / 2, 15}, {ll / 2, -4 + ll / 2, 16}, {ll / 2, -4 + ll / 2, 17}, {ll / 2, -4 + ll / 2, 18}, {ll / 2, -4 + ll / 2, 19}, {ll / 2, -3 + ll / 2, 11}, {ll / 2, -3 + ll / 2, 12}, {ll / 2, -3 + ll / 2, 13}, {ll / 2, -3 + ll / 2, 14}, {ll / 2, -3 + ll / 2, 15}, {ll / 2, -3 + ll / 2, 16}, {ll / 2, -3 + ll / 2, 17}, {ll / 2, -3 + ll / 2, 18}, {ll / 2, -3 + ll / 2, 19}, {ll / 2, -2 + ll / 2, 11}, {ll / 2, -2 + ll / 2, 12}, {ll / 2, -2 + ll / 2, 13}, {ll / 2, -2 + ll / 2, 14}, {ll / 2, -2 + ll / 2, 15}, {ll / 2, -2 + ll / 2, 16}, {ll / 2, -2 + ll / 2, 17}, {ll / 2, -2 + ll / 2, 18}, {ll / 2, -2 + ll / 2, 19}, {ll / 2, -1 + ll / 2, 11}, {ll / 2, -1 + ll / 2, 12}, {ll / 2, -1 + ll / 2, 13}, {ll / 2, -1 + ll / 2, 14}, {ll / 2, -1 + ll / 2, 15}, {ll / 2, -1 + ll / 2, 16}, {ll / 2, -1 + ll / 2, 17}, {ll / 2, -1 + ll / 2, 18}, {ll / 2, -1 + ll / 2, 19}, {ll / 2, 1 + ll / 2, 11}, {ll / 2, 1 + ll / 2, 12}, {ll / 2, 1 + ll / 2, 13}, {ll / 2, 1 + ll / 2, 14}, {ll / 2, 1 + ll / 2, 15}, {ll / 2, 1 + ll / 2, 16}, {ll / 2, 1 + ll / 2, 17}, {ll / 2, 1 + ll / 2, 18}, {ll / 2, 1 + ll / 2, 19}, {ll / 2, 2 + ll / 2, 11}, {ll / 2, 2 + ll / 2, 12}, {ll / 2, 2 + ll / 2, 13}, {ll / 2, 2 + ll / 2, 14}, {ll / 2, 2 + ll / 2, 15}, {ll / 2, 2 + ll / 2, 16}, {ll / 2, 2 + ll / 2, 17}, {ll / 2, 2 + ll / 2, 18}, {ll / 2, 2 + ll / 2, 19}, {ll / 2, 3 + ll / 2, 11}, {ll / 2, 3 + ll / 2, 12}, {ll / 2, 3 + ll / 2, 13}, {ll / 2, 3 + ll / 2, 14}, {ll / 2, 3 + ll / 2, 15}, {ll / 2, 3 + ll / 2, 16}, {ll / 2, 3 + ll / 2, 17}, {ll / 2, 3 + ll / 2, 18}, {ll / 2, 3 + ll / 2, 19}, {ll / 2, 4 + ll / 2, 11}, {ll / 2, 4 + ll / 2, 12}, {ll / 2, 4 + ll / 2, 13}, {ll / 2, 4 + ll / 2, 14}, {ll / 2, 4 + ll / 2, 15}, {ll / 2, 4 + ll / 2, 16}, {ll / 2, 4 + ll / 2, 17}, {ll / 2, 4 + ll / 2, 18}, {ll / 2, 4 + ll / 2, 19}, {ll / 2, 5 + ll / 2, 11}, {ll / 2, 5 + ll / 2, 12}, {ll / 2, 5 + ll / 2, 13}, {ll / 2, 5 + ll / 2, 14}, {ll / 2, 5 + ll / 2, 15}, {ll / 2, 5 + ll / 2, 16}, {ll / 2, 5 + ll / 2, 17}, {ll / 2, 5 + ll / 2, 18}, {ll / 2, 5 + ll / 2, 19}, {-1 + ll / 2, ll / 2, 1}, {-1 + ll / 2, ll / 2, 2}, {-1 + ll / 2, ll / 2, 3}, {-1 + ll / 2, ll / 2, 4}, {-1 + ll / 2, ll / 2, 5}, {-1 + ll / 2, ll / 2, 6}, {-1 + ll / 2, ll / 2, 7}, {-1 + ll / 2, ll / 2, 8}, {-1 + ll / 2, ll / 2, 9}, {-1 + ll / 2, ll / 2, 10}, {1 + ll / 2, ll / 2, 1}, {1 + ll / 2, ll / 2, 2}, {1 + ll / 2, ll / 2, 3}, {1 + ll / 2, ll / 2, 4}, {1 + ll / 2, ll / 2, 5}, {1 + ll / 2, ll / 2, 6}, {1 + ll / 2, ll / 2, 7}, {1 + ll / 2, ll / 2, 8}, {1 + ll / 2, ll / 2, 9}, {1 + ll / 2, ll / 2, 10}, {ll / 2, -1 + ll / 2, 1}, {ll / 2, -1 + ll / 2, 2}, {ll / 2, -1 + ll / 2, 3}, {ll / 2, -1 + ll / 2, 4}, {ll / 2, -1 + ll / 2, 5}, {ll / 2, -1 + ll / 2, 6}, {ll / 2, -1 + ll / 2, 7}, {ll / 2, -1 + ll / 2, 8}, {ll / 2, -1 + ll / 2, 9}, {ll / 2, -1 + ll / 2, 10}, {ll / 2, 1 + ll / 2, 1}, {ll / 2, 1 + ll / 2, 2}, {ll / 2, 1 + ll / 2, 3}, {ll / 2, 1 + ll / 2, 4}, {ll / 2, 1 + ll / 2, 5}, {ll / 2, 1 + ll / 2, 6}, {ll / 2, 1 + ll / 2, 7}, {ll / 2, 1 + ll / 2, 8}, {ll / 2, 1 + ll / 2, 9}, {ll / 2, 1 + ll / 2, 10}}};
    
        matrix<double> newdat(Nfan, 3);
    for (int i = 0; i < Nfan; i++)
    {
        newdat(i, 0) = fan_data[i][0];
        newdat(i, 1) = fan_data[i][1];
        newdat(i, 2) = fan_data[i][2];
    }


    matrix<int> fanpairs(388,2);
    std::array<std::array<int, 2>, 388> fan_index = {{{0, 1}, {0, 9}, {1, 2}, {1, 10}, {2, 3}, {2, 11}, {3, 4}, {3, 12}, {4, 5}, {4, 13}, {5, 6}, {5, 14}, {6, 7}, {6, 15}, {7, 8}, {7, 16}, {8, 17}, {9, 10}, {9, 18}, {10, 11}, {10, 19}, {11, 12}, {11, 20}, {12, 13}, {12, 21}, {13, 14}, {13, 22}, {14, 15}, {14, 23}, {15, 16}, {15, 24}, {16, 17}, {16, 25}, {17, 26}, {18, 19}, {18, 27}, {19, 20}, {19, 28}, {20, 21}, {20, 29}, {21, 22}, {21, 30}, {22, 23}, {22, 31}, {23, 24}, {23, 32}, {24, 25}, {24, 33}, {25, 26}, {25, 34}, {26, 35}, {27, 28}, {27, 36}, {28, 29}, {28, 37}, {29, 30}, {29, 38}, {30, 31}, {30, 39}, {31, 32}, {31, 40}, {32, 33}, {32, 41}, {33, 34}, {33, 42}, {34, 35}, {34, 43}, {35, 44}, {36, 37}, {36, 45}, {36, 198}, {37, 38}, {37, 46}, {38, 39}, {38, 47}, {39, 40}, {39, 48}, {40, 41}, {40, 49}, {41, 42}, {41, 50}, {42, 43}, {42, 51}, {43, 44}, {43, 52}, {44, 53}, {45, 46}, {45, 54}, {45, 135}, {45, 144}, {46, 47}, {46, 55}, {46, 136}, {46, 145}, {47, 48}, {47, 56}, {47, 137}, {47, 146}, {48, 49}, {48, 57}, {48, 138}, {48, 147}, {49, 50}, {49, 58}, {49, 139}, {49, 148}, {50, 51}, {50, 59}, {50, 140}, {50, 149}, {51, 52}, {51, 60}, {51, 141}, {51, 150}, {52, 53}, {52, 61}, {52, 142}, {52, 151}, {53, 62}, {53, 143}, {53, 152}, {54, 55}, {54, 63}, {54, 208}, {55, 56}, {55, 64}, {56, 57}, {56, 65}, {57, 58}, {57, 66}, {58, 59}, {58, 67}, {59, 60}, {59, 68}, {60, 61}, {60, 69}, {61, 62}, {61, 70}, {62, 71}, {63, 64}, {63, 72}, {64, 65}, {64, 73}, {65, 66}, {65, 74}, {66, 67}, {66, 75}, {67, 68}, {67, 76}, {68, 69}, {68, 77}, {69, 70}, {69, 78}, {70, 71}, {70, 79}, {71, 80}, {72, 73}, {72, 81}, {73, 74}, {73, 82}, {74, 75}, {74, 83}, {75, 76}, {75, 84}, {76, 77}, {76, 85}, {77, 78}, {77, 86}, {78, 79}, {78, 87}, {79, 80}, {79, 88}, {80, 89}, {81, 82}, {81, 90}, {82, 83}, {82, 91}, {83, 84}, {83, 92}, {84, 85}, {84, 93}, {85, 86}, {85, 94}, {86, 87}, {86, 95}, {87, 88}, {87, 96}, {88, 89}, {88, 97}, {89, 98}, {90, 91}, {91, 92}, {92, 93}, {93, 94}, {94, 95}, {95, 96}, {96, 97}, {97, 98}, {99, 100}, {99, 108}, {100, 101}, {100, 109}, {101, 102}, {101, 110}, {102, 103}, {102, 111}, {103, 104}, {103, 112}, {104, 105}, {104, 113}, {105, 106}, {105, 114}, {106, 107}, {106, 115}, {107, 116}, {108, 109}, {108, 117}, {109, 110}, {109, 118}, {110, 111}, {110, 119}, {111, 112}, {111, 120}, {112, 113}, {112, 121}, {113, 114}, {113, 122}, {114, 115}, {114, 123}, {115, 116}, {115, 124}, {116, 125}, {117, 118}, {117, 126}, {118, 119}, {118, 127}, {119, 120}, {119, 128}, {120, 121}, {120, 129}, {121, 122}, {121, 130}, {122, 123}, {122, 131}, {123, 124}, {123, 132}, {124, 125}, {124, 133}, {125, 134}, {126, 127}, {126, 135}, {127, 128}, {127, 136}, {128, 129}, {128, 137}, {129, 130}, {129, 138}, {130, 131}, {130, 139}, {131, 132}, {131, 140}, {132, 133}, {132, 141}, {133, 134}, {133, 142}, {134, 143}, {135, 136}, {135, 218}, {136, 137}, {137, 138}, {138, 139}, {139, 140}, {140, 141}, {141, 142}, {142, 143}, {144, 145}, {144, 153}, {144, 228}, {145, 146}, {145, 154}, {146, 147}, {146, 155}, {147, 148}, {147, 156}, {148, 149}, {148, 157}, {149, 150}, {149, 158}, {150, 151}, {150, 159}, {151, 152}, {151, 160}, {152, 161}, {153, 154}, {153, 162}, {154, 155}, {154, 163}, {155, 156}, {155, 164}, {156, 157}, {156, 165}, {157, 158}, {157, 166}, {158, 159}, {158, 167}, {159, 160}, {159, 168}, {160, 161}, {160, 169}, {161, 170}, {162, 163}, {162, 171}, {163, 164}, {163, 172}, {164, 165}, {164, 173}, {165, 166}, {165, 174}, {166, 167}, {166, 175}, {167, 168}, {167, 176}, {168, 169}, {168, 177}, {169, 170}, {169, 178}, {170, 179}, {171, 172}, {171, 180}, {172, 173}, {172, 181}, {173, 174}, {173, 182}, {174, 175}, {174, 183}, {175, 176}, {175, 184}, {176, 177}, {176, 185}, {177, 178}, {177, 186}, {178, 179}, {178, 187}, {179, 188}, {180, 181}, {181, 182}, {182, 183}, {183, 184}, {184, 185}, {185, 186}, {186, 187}, {187, 188}, {189, 190}, {190, 191}, {191, 192}, {192, 193}, {193, 194}, {194, 195}, {195, 196}, {196, 197}, {197, 198}, {199, 200}, {200, 201}, {201, 202}, {202, 203}, {203, 204}, {204, 205}, {205, 206}, {206, 207}, {207, 208}, {209, 210}, {210, 211}, {211, 212}, {212, 213}, {213, 214}, {214, 215}, {215, 216}, {216, 217}, {217, 218}, {219, 220}, {220, 221}, {221, 222}, {222, 223}, {223, 224}, {224, 225}, {225, 226}, {226, 227}, {227, 228}}};
    
    for (int i = 0; i < 388; i++)
    {
        (fanpairs)(i, 0) = fan_index[i][0];
        (fanpairs)(i, 1) = fan_index[i][1];
    }

    ShellProperties B;
    B.k = 100.;
    B.rm = 0.;
    B.par = fanpairs;
    B.posi =  newdat;

    anchoringpointsharmonic a1;
    anchoringpointsharmonic a2;
    anchoringpointsharmonic a3;
    anchoringpointsharmonic a4;
    a1.epsilon = 100.;
    a2.epsilon = 100.;
    a3.epsilon = 100.;
    a4.epsilon = 100.;
    a1.anchorpointx = -2 + ll / 2.;
    a1.anchorpointy = ll / 2;
    a1.anchorpointz = 0;

    a2.anchorpointx = 2 + ll / 2.;
    a2.anchorpointy = ll / 2;
    a2.anchorpointz = 0;

    a3.anchorpointx = ll / 2.;
    a3.anchorpointy = -2 + ll / 2;
    a3.anchorpointz = 0;

    a4.anchorpointx = ll / 2.;
    a4.anchorpointy = 2 + ll / 2;
    a4.anchorpointz = 0;
    vector<anchoringpointsharmonic> ap;
    ap.push_back(a1);
    ap.push_back(a2);
    ap.push_back(a3);
    ap.push_back(a4);

    vector<int> index_anchored;
    index_anchored.push_back(189);
    index_anchored.push_back(199);
    index_anchored.push_back(209);
    index_anchored.push_back(219);

    B.anchors = ap;
    B.index_anchored = index_anchored;



    bool cc;
    NanotubeAssembly A(ll, (g.no_types)*5000+1,cc);



    double interactionenergy = 100.;


    GeneralPatch c(CreateGeneralPatchFan(interactionenergy, 1, 1.2, 0.6, g,Nfan));


    A.setpots(c);
    A.setkT(1.0);
    A.setviscosity(1.0);
    A.run_fan_equil(20000000, 1000, g, B, "");

    
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