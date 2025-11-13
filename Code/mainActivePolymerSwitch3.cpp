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
#include <chrono>
#include <iomanip>
// #include <sys/ioctl.h>
// #include <fcntl.h>
// #include <time.h>
// #include <sys/time.h>
#include <sys/stat.h>
#include <random>
#include <algorithm>
#include <parallel/algorithm>
// #include <thrust/host_vector.h>
// #include <thrust/device_vector.h>
#include <unistd.h>
#include <execinfo.h>
#include <signal.h>

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
#include "Condensate/Condensate.h"
#include "NCGas/NCGas.h"

// #include "NCGasR.h"
// #include "Microtubule.h"

// #include "MDGPU.cu"

using namespace std;

matrix<double> harmoniccurl(matrix<double> &dat, double m)
{
    matrix<double> F2(dat.getnrows(), dat.getncols());

    for (int i = 0; i < dat.getnrows() - 1; i++)
    {
        double x0 = dat(i, 0);
        double y0 = dat(i, 1);
        double z0 = dat(i, 2);

        double x1 = dat(i + 1, 0);
        double y1 = dat(i + 1, 1);
        double z1 = dat(i + 1, 2);

        F2(i + 1, 0) += m * (-z0 + z1);
        F2(i + 1, 1) += m * (-x0 + x1);
        F2(i + 1, 2) += m * (-y0 + y1);

        F2(i, 0) += -m * (-z0 + z1);
        F2(i, 1) += -m * (-x0 + x1);
        F2(i, 2) += -m * (-y0 + y1);
    }
    return F2;
}

int main(int argc, char **argv)
{
    unsigned long seed = mix(clock(), time(NULL), getpid());
    srand(seed);

    // srand(time(NULL));

    // uint64_t microseconds_since_epoch = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // // cout << microseconds_since_epoch << endl;
    // // cout << time(NULL) << endl;
    // int seed = microseconds_since_epoch % time(NULL);
    // srand(seed);
    double var = 1.;
    // double var;
    // int num;
    // if (argc == 3)
    // {
    //     var = atof(argv[1]);
    //     num = atof(argv[2]);
    // }
    // else
    // {
    //     error("not enough arguments");
    // }

    int dimension = 3;

    double ll = 100.0;
    vector1<bool> is_periodic(dimension, true);
    cube bc(ll, is_periodic, dimension);

    int ccc;
    int num2;

    double att_epp = 1.0;

    // HarmonicPotential pot(att_epp, 0.);

    FENEPotential pot(30.0, 1.5);
    LJPotential pot2(att_epp, 1.0);

    // double att_epp2 = 0.1;

    // QuarticPotential pot2(att_epp2, 0.);

    LangevinNVT b(bc);

    double kT = 1.0;
    double dt = 0.005;
    double eta = 20.;
    // gamma = eta;

    // b.setinteractions(nswca);
    b.setkT(kT);
    b.setdt(dt);
    b.setgamma(eta);
    b.setm(1.0);

    int NN = 1000;

    double baseline = (ll / 2.);

    matrix<double> dat(NN, dimension);
    matrix<double> moms(NN, dimension);

    // place polymer, we can use the grid
    string filename = "final.csv";
    double T;
    bool err;
    matrix<double> newpos = importcsv(filename, T, err);

    //

    dat = newpos;

    b.setdat(dat);
    b.setmom(moms);

    b.setgeometry(bc);

    num2 = (int)(100. / 4);
    matrix<int> boxes = bc.generate_boxes_relationships(num2, ccc);

    int runtime = 20000000;
    // double tem[10] = {0.44753, 0.150461, 0.477014, 0.281367, -0.300378, 0.352483, -0.277297, -0.427009, 0.144324, 0.212938};
    vector1<double> tempmod(NN, 1.);

    // vector1<double> createrandom(NN);
    // for (int i = 0; i < NN; i++)
    // {
    //     createrandom[i] = (2. * (double)rand() / (double)(RAND_MAX)-1.);
    // }

    // double meani = meanish(createrandom);

    // for(int i = 0  ; i < NN ; i++ ) {
    //     createrandom[i] -= meani;
    // }

    matrix<int> bp2(NN - 1, 2);
    for (int i = 0; i < NN - 1; i++)
    {
        bp2(i, 0) = i;
        bp2(i, 1) = i + 1;
    }

    // double diff[9] = {5.52701, 9.84486, 1.11208, 2.52656, 2.06907, 2.405, 2.6374, 6.70836, 6.73533};

    vector1<double> mags(NN - 1);
    for (int i = 0; i < NN - 1; i++)
    {
        mags[i] = 1.;
    }

    matrix<double> F(NN, dimension);
    // matrix<double> F2(NN, dimension);
    double bounding_pot = 0.0;
    int every = 1000;

    stringstream ss1;
    stringstream ss2;

    ss1 << var;
    ss2 << NN;

    string vs = ss1.str();
    string ns = ss2.str();

    string ext = ".csv";

    string temps = "temps_var=" + vs + "_N=" + ns + ext;
    string datas = "data_var=" + vs + "_N=" + ns + ext;

    ofstream mytemps;
    mytemps.open(temps.c_str());
    // mytemps <<= tempmod;
    // mytemps.close();

    ofstream myfile;
    myfile.open(datas.c_str());

    for (int j = 450; j < 550; j++)
    {

        tempmod[j] = 1.;
    }

    // from 0 to 600 is excitable
    

    double rate_on = 0.001 * 200.;
    double rate_off = 0.0001 * 200;

    int block_size = 6;
    int potential_blocks = floor(NN/block_size);
    vector<int> switchable(potential_blocks);
    for(int i = 0  ; i < potential_blocks  ; i++) {
        switchable[i] = i;
    }

    double proportion_switchable = 0.6;
    int tot_monomers_switchable =  0.6*NN;
    int block_switchable = tot_monomers_switchable/block_size;

    std::vector<int> out;
    

    std::random_device rd;
    std::mt19937 gen(rd());

    std::sample(switchable.begin(), switchable.end(), std::back_inserter(out), block_switchable, gen);

    matrix<int> pos(block_switchable,block_size);

    for(int i = 0  ; i < block_switchable ; i++) {
        int j = out[i]*block_size;

        for(int k  = 0 ; k < block_size ; k++) {
            pos(i,k) =  j+k;
        }
    }

    vector1<bool> excitation_on(block_switchable, false);



    
    int blocks_on = 10;


    matrix<int> *pairs = b.calculatepairs_parallel(boxes, 3.5);

    for (int i = 0; i < runtime; i++)
    {
        if (i % 200 == 0)
        {

            double r = (double)rand() / (double)RAND_MAX;
            // int len = rand() % 60;
            int reg = rand() % block_switchable;
            int alls = total_bool(excitation_on);
            if (r < rate_on && excitation_on[reg] == false && alls <= blocks_on)
            {

                for (int pp = 0; pp < block_size; pp++)
                {
                    int mon_iter = pos(reg,pp);
                    tempmod[mon_iter] = 2.5;
                }
                excitation_on[reg] = true;
            }
            r = (double)rand() / (double)RAND_MAX;
            if (r < rate_off)
            {
                vector<int> pos_on;
                for (int j = 0; j < block_switchable; j++)
                {
                    if (excitation_on[j] == true)
                        pos_on.push_back(j);
                }
                if (pos_on.size() > 0)
                {
                    int reg = pos[rand() % pos_on.size()];
                    for (int pp = 0; pp < block_size; pp++)
                    {
                        int mon_iter = pos(reg, pp);
                        tempmod[mon_iter] = 1.0;
                    }
                    excitation_on[reg] = false;
                }
            }
        }

        // if(i%1000000 == 0 && i > 0 ) {

        //         string stra = "temps";
        //         stringstream ssx;
        //         ssx << (i/1000000);
        //         stra += ssx.str();
        //         outfunc(tempmod, stra);

        // }

        cout << "i:" << i << endl;
        if (i > 0 && i % 20 == 0)
        {
            delete pairs;
            pairs = b.calculatepairs_parallel(boxes, 3.5);
        }
        // cout << i << endl;
        matrix<double> F1(b.calculateforces(bp2, pot));

        matrix<double> F2(b.calculateforces(*pairs, pot2));

        // matrix<double> F2(b.calculateforcesDV(bp2, pot2, mags));

        // for(int i1 = 0 ; i1 < NN ; i1++) {
        //     for(int j1 = 0  ; j1 < dimension ; j1++) {
        //         F2(i1,j1) = -bounding_pot*(b.getcoordinate(i1,j1)-ll/2.);
        //     }
        // }
        // matrix<double> pos2 = b.getdat();
        // matrix<double> F2 = harmoniccurl(pos2, mf);
        F = F1; // + F2;
        F += F2;

        matrix<double> R(NN, dimension);

        for (int i1 = 0; i1 < NN; i1++)
        {
            for (int j = 0; j < dimension; j++)
            {
                R(i1, j) = tempmod[i1] * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
            }
        }
        if (i % every == 0)
        {
            cout << i << endl;
            matrix<double> pos = b.getdat();
            myfile <<= pos;
            mytemps <<= tempmod;
            mytemps << endl;
        }
        b.advance_mom(F, R);

        b.advance_pos();
    }
    myfile.close();

    return 0;
}