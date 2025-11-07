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
    opt1[0]=1;
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

    vector1<int> opt5(6,1);
    opt5[5]=0;

    vector1<int> opt6(6,1);

    // cout << opt1 << endl;


    std::vector<int> weights = {1, 4, 3, 2, 1, 1};

    // random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // discrete distribution with given weights
    std::discrete_distribution<> dist(weights.begin(), weights.end());

    matrix<int> olg(str1,7);
    for(int i = 0  ; i < str1 ; i++) {
        
        int n = dist(gen) + 1;
        olg(i, 0) = n; // between 1 and 6


        //std::vector<int> v = {1, 2, 3, 4, 5, 6};

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


        switch(n) {
            case 1 : {
                for(int j = 0 ; j < 6 ; j++)
                    olg(i,j+1) = opt1[j];
 
                break;
                }
            case 2 : {
                int randr = rand() % 2;
                if(randr ==0) {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt2_1[j];
                }
                else{
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt2_2[j];
                }

                break;
            }
            case 3 : {
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
            case 4 : {
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
            case 5 : {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt5[j];

                break;
            }
            case 6 : {
                for (int j = 0; j < 6; j++)
                    olg(i, j + 1) = opt6[j];

                break;
            }
        }


        string result;

        result += std::to_string(n);
        for(int j =1 ; j < 7 ; j++) {
            result += std::to_string(olg(i,j));
        }
        s.push_back(result);
    }


    for(int i = 0 ; i < str1 ; i++) {
        for(int j = i ; j < str1 ; j++) {
            int totn = (olg(i, 0) * olg(j, 0));
            vector1<int> ints(totn,1);
            for(int k  = 0 ; k < totn ; k++) {
                int ran = rand() % 6;
                if(ran==5) ints[k] = 0;
            }
            string result;

            for(int k = 0 ; k < totn ; k++) {
                result += to_string(ints[k]);
            }
            s.push_back(result);

        }
    }

    s.push_back("-5");

    vector1<int> weig(str1);

    for(int i = 0 ; i < str1 ; i++) {
        weig[i]= 1 + rand() % 3;
    }

    string ress;

    for(int i = 0  ; i < str1 ; i++) {
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
    for(int i = 0  ; i < s.size() ; i++) {
        myfile << s[i] << endl;
    }
    myfile.close();

    
    


    geneticcode g(s);

    // g.rate = 0.0;



    bool cc;
    NanotubeAssembly A(20., (g.no_types)*5000+1,cc);
    

    // cout << *(g.patch_num) << endl;
    // for(int i = 0  ; i < g.interactions.size() ; i++) {
    // for(int k = 0 ; k < g.interactions[i].size() ; k++ ) {
    // cout << g.interactions[i][k] << " ";
    // }
    // cout << endl;
    // }
    // pausel();

    GeneralPatch c(CreateGeneralPatch(100., 1, 1.2, 0.6, g));


    
    A.setpots(c);
    A.setkT(1.0);
    A.setviscosity(1.0);

    A.run_box_equil(10000000,1000,100.,g,"");
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