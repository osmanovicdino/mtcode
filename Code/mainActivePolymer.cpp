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

    // srand(time(NULL));

    uint64_t microseconds_since_epoch = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // cout << microseconds_since_epoch << endl;
    // cout << time(NULL) << endl;
    int seed = microseconds_since_epoch % time(NULL);
    srand(seed);
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

    FENEPotential pot(30.0,1.5);
    LJPotential pot2(att_epp,1.0);

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

    dat(0,0) = ll/2.;
    dat(0,1) = ll/2.;
    dat(0,2) = ll/2;

    //place polymer, we can use the grid

    vector<double> possible_pos_x;
    vector<double> possible_pos_y;
    vector<double> possible_pos_z;


    int pp = 100;
    possible_pos_x.reserve(pp);
    possible_pos_y.reserve(pp);
    possible_pos_z.reserve(pp);

    for (int j = 0; j < pp; j++)
    {
        double x = 0.5 +  j;
        // vector<double> b;
        possible_pos_x.push_back(x);
        possible_pos_y.push_back(x);
        possible_pos_z.push_back(x);
        // possible_pos.push_back(b);
    }

    struct tripletstorage {
        int i;
        int j;
        int k;
        tripletstorage(int ii,int jj,int kk) {
            i =  ii;
            j =  jj;
            k =  kk;
        }
        tripletstorage(const tripletstorage &c) {
            i = c.i;
            j = c.j;
            k = c.k;
        }
        tripletstorage& operator=(const tripletstorage &c) {
            i = c.i;
            j = c.j;
            k = c.k;

            return *this;
        }

        bool operator==(const tripletstorage &b) {
            if (i == b.i && j == b.j && k == b.k) {
                return true;
            }
            else return false;
        }
    };

    int si = 50;// rand() % 100;
    int sj = 50; //rand() % 100;
    int sk = 50; //rand() % 100;

    vector<tripletstorage> pos_iter;
    matrix<double> newpos(NN,3);

    newpos(0, 0) = possible_pos_x[si];
    newpos(0, 1) = possible_pos_x[sj];
    newpos(0, 2) = possible_pos_x[sk];

    tripletstorage temp(si, sj, sk);
    pos_iter.push_back(temp);

    tripletstorage c1(-1,0,0);
    tripletstorage c2(1, 0, 0);
    tripletstorage c3(0, -1, 0);
    tripletstorage c4(0, 1, 0);
    tripletstorage c5(0, 0, -1);
    tripletstorage c6(0, 0, 1);

    vector<tripletstorage> ts;
    ts.push_back(c1);
    ts.push_back(c2);
    ts.push_back(c3);
    ts.push_back(c4);
    ts.push_back(c5);
    ts.push_back(c6);
    int sxn = si;
    int syn = sj;
    int szn = sk;
    for (int i = 1; i < NN; i++)
    {
       // cout << i << endl;
        bool overlap = true;
        tripletstorage temp2(0, 0, 0);

        while(overlap) {
        int direction = rand() % 6;
        tripletstorage mot = ts[direction];
        


        int sxnn = sxn + mot.i;
        int synn = syn + mot.j;
        int sznn = szn + mot.k;



        if(sxnn < 0) sxnn = 99;
        if(sxnn > 99) sxnn = 0;


        if(synn < 0) synn = 99;
        if(synn > 99) synn  = 0;


        if(sznn < 0) sznn = 99;
        if(sznn > 99) sznn  = 0;

        temp2 = tripletstorage(sxnn, synn, sznn);

        for(int j = 0  ; j < pos_iter.size() ; j++) {
            if(temp2 ==  pos_iter[j] ) { overlap = true; break; }
            
            else {
                overlap = false;

            }

            }
        }

        sxn = temp2.i;
        syn = temp2.j;
        szn = temp2.k;

        pos_iter.push_back(temp2);
        newpos(i, 0) = possible_pos_x[sxn];
        newpos(i, 1) = possible_pos_x[syn];
        newpos(i, 2) = possible_pos_x[szn];
    }


    //

    dat = newpos;

    b.setdat(dat);
    b.setmom(moms);

    b.setgeometry(bc);

    num2 = (int)(100./4);
    matrix<int> boxes = bc.generate_boxes_relationships(num2, ccc);

    int runtime = 10000000;
    // double tem[10] = {0.44753, 0.150461, 0.477014, 0.281367, -0.300378, 0.352483, -0.277297, -0.427009, 0.144324, 0.212938};
    vector1<double> tempmod(NN,1.);

    // vector1<double> createrandom(NN);
    // for (int i = 0; i < NN; i++)
    // {
    //     createrandom[i] = (2. * (double)rand() / (double)(RAND_MAX)-1.);
    // }

    // double meani = meanish(createrandom);

    // for(int i = 0  ; i < NN ; i++ ) {
    //     createrandom[i] -= meani;
    // }

    for (int i = 0; i < 100; i++)
    {

        tempmod[i] = 1.;
    }

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
    //matrix<double> F2(NN, dimension);
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
    mytemps <<= tempmod;
    mytemps.close();

    ofstream myfile;
    myfile.open(datas.c_str());




    matrix<int> *pairs = b.calculatepairs_parallel(boxes, 3.5);

    for (int i = 0; i < runtime; i++)
    {
        cout << "i:" << i << endl;
        if (i > 0 && i % 20 == 0)
        {
            delete pairs;
            pairs = b.calculatepairs_parallel(boxes, 3.5);
        }
        // cout << i << endl;
        matrix<double> F1(b.calculateforces(bp2, pot));



        matrix<double> F2(b.calculateforces(*pairs,pot2));

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
        }
        b.advance_mom(F, R);

        b.advance_pos();
    }
    myfile.close();

    return 0;
}