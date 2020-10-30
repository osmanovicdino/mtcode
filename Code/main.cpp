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
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "Basic/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"
#include "MDBase/potential.h"
#include "MDBase/MD.h"
#include "MDBase/Langevin.h"
#include "MDBase/LangevinR.h"
#include "Condensate/Condensate.h"

// #include "NCGasR.h"
// #include "Microtubule.h"



//#include "MDGPU.cu"

using namespace std;



int main(int argc, char** argv) {

srand (time(NULL));

Condensate A(25.0, 1000, 2);


vector1<double> as(16,10.0);
vector1<double> bs(16,2.);
vector1<double> cs(16,pi/3.);


//as[0] = 100.0;



set_potential_bundle_tetrahedral(A, as, bs, cs);

for(double kT = 1.0 ; kT > 0.49 ; kT -= 0.1) {

A.obj->setkT(kT);
stringstream ss;
ss << kT;
    
string base = "_kT=";
base += ss.str();


A.run(1000000,1000,base);
}



// cout << (A.obj)->getdat() << endl;

// cout << (A.obj)->getorientation() << endl;

/*
vector1<bool> pb(3,true);
cube geo(5.1,pb,3);

LangevinNVTR a(geo);


matrix<double> dat(10,3);

//for the given system we can rotate everything

//matrix<double> Rot(3,3);

// double thetax  =pi*(double)rand()/(double)RAND_MAX;
 
// double thetay = pi * (double)rand() / (double)RAND_MAX;

// double thetaz = pi * (double)rand() / (double)RAND_MAX;


// double thetax = 0.0;

// double thetay = 0.0;

// double thetaz = 0.0;

// Rot(0, 0) = cos(thetay) * cos(thetaz);
// Rot(0, 1) = -(cos(thetay) * sin(thetaz));
// Rot(0, 2) = sin(thetay);
// Rot(1, 0) = cos(thetaz) * sin(thetax) * sin(thetay) + cos(thetax) * sin(thetaz);
// Rot(1, 1) = cos(thetax) * cos(thetaz) - sin(thetax) * sin(thetay) * sin(thetaz);
// Rot(1, 2) = -(cos(thetay) * sin(thetax));
// Rot(2, 0) = -(cos(thetax) * cos(thetaz) * sin(thetay)) + sin(thetax) * sin(thetaz);
// Rot(2, 1) = cos(thetaz) * sin(thetax) + cos(thetax) * sin(thetay) * sin(thetaz);
// Rot(2, 2) = cos(thetax) * cos(thetay);

// double nx = 1;
// double ny = 0;
// double nz = 0;
// double theta = -0.2;
// double phi = pi/2+.2;

// double nxb = cos(theta) * sin(phi);
// double nyb = sin(theta) * sin(phi);
// double nzb = cos(phi);

// double nx = Rot(0, 0) * nxb + Rot(0, 1) * nyb + Rot(0,2) * nzb;

// double ny = Rot(1, 0) * nxb + Rot(1, 1) * nyb + Rot(1, 2) * nzb;

// double nz = Rot(2, 0) * nxb + Rot(2, 1) * nyb + Rot(2, 2) * nzb;

double nx = 1;
double ny = 0;
double nz = 0;

KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(nx, ny, nz, nx, ny, nz, 100.0, 2., pi / 3., 0.75);
KernFrenkelOnePatch2 *pot2 = new KernFrenkelOnePatch2(nx, ny, nz, -nx, ny, nz, 100.0, 2., pi / 3., 0.75);
KernFrenkelOnePatch2 *pot3 = new KernFrenkelOnePatch2(-nx, ny, nz, nx, ny, nz, 100.0, 2., pi / 3., 0.75);
KernFrenkelOnePatch2 *pot4 = new KernFrenkelOnePatch2(-nx, ny, nz, -nx, ny, nz, 100.0, 2., pi / 3., 0.75);

cout << "there" << endl;

vector1<potentialtheta3D*> pots(4);

cout << "here" << endl;
pots[0] = pot1;
pots[1] = pot2;
pots[2] = pot3;
pots[3] = pot4;

cout << "here" << endl;

int tot = 25;

matrix<double> pos(tot,3);

int iter = 0;
for(int i = 0  ; i < 5 ; i++ ) {
    for (int j = 0; j < 5; j++)
    {
        //for(int j = 0  ; j < 3 ; j++) {
        pos(iter,0) =i+0.05;
        pos(iter,1) =j+0.05;
        pos(iter,2) = 1.;
        iter++;
    }
}




// double pos1x = Rot(0, 0) * pos(0, 0) + Rot(0, 1) * pos(0, 1) + Rot(0, 2) * pos(0, 2);

// double pos1y = Rot(1, 0) * pos(0, 0) + Rot(1, 1) * pos(0, 1) + Rot(1, 2) * pos(0, 2);

// double pos1z = Rot(2, 0) * pos(0, 0) + Rot(2, 1) * pos(0, 1) + Rot(2, 2) * pos(0, 2);

// double pos2x = Rot(0, 0) * pos(1, 0) + Rot(0, 1) * pos(1, 1) + Rot(0, 2) * pos(1, 2);

// double pos2y = Rot(1, 0) * pos(1, 0) + Rot(1, 1) * pos(1, 1) + Rot(1, 2) * pos(1, 2);

// double pos2z = Rot(2, 0) * pos(1, 0) + Rot(2, 1) * pos(1, 1) + Rot(2, 2) * pos(1, 2);

// double tx = (pos1x - 50.0);

// double ty = (pos1y - 50.0);

// double tz = (pos1z - 50.0);

// pos(0, 0) = pos1x - tx;
// pos(0, 1) = pos1y - ty;
// pos(0, 2) = pos1z - tz;
// pos(1, 0) = pos2x - tx;
// pos(1, 1) = pos2y - ty;
// pos(1, 2) = pos2z - tz;

// cout << Rot << endl;
// cout << pos << endl;
// pausel();

// matrix<double> mom(2,3);
// matrix<double> amo(2,3);

// matrix<double> ori(2,9);

// ori(0, 0) = 1;
// ori(0, 1) = 0;
// ori(0, 2) = 0;
// ori(0, 3) = 0;
// ori(0, 4) = 1;
// ori(0, 5) = 0;
// ori(0, 6) = 0;
// ori(0, 7) = 0;
// ori(0, 8) = 1;

// ori(1, 0) = 1;
// ori(1, 1) = 0;
// ori(1, 2) = 0;
// ori(1, 3) = 0;
// ori(1, 4) = 1;
// ori(1, 5) = 0;
// ori(1, 6) = 0;
// ori(1, 7) = 0;
// ori(1, 8) = 1;


a.initialize(pos);

matrix<double> F(tot, 3);
matrix<double> T(tot, 3);
matrix<double> zeromatrix(tot,3);
// matrix<int> pairs(1,2);
// pairs(0,0) = 0;
// pairs(0,1) = 1;


WCAPotential wsa(1.0,1.0,0.0);


a.setdt(0.001);

double viscosity = 10.0;
double hdradius = 0.5;
double fac1 = 1.;
double fac2 = 1.;
a.setgamma(fac1*6*pi*viscosity*hdradius);
a.setgammar(fac2*8*pi*viscosity*hdradius*hdradius*hdradius);
a.setkT(1.0);

int ccc;
int num = 2;
matrix<int> boxes = (a).getgeo().generate_boxes_relationships(num, ccc);

matrix<int> *pairs = a.calculatepairs(boxes,3.5);

a.calculateforces(*pairs, wsa);
a.calculate_forces_and_torques3D(*pairs, pots, F, T);

// matrix<double> tempT = T;


a.create_forces_and_torques_sphere(F, T);



// ofstream myfile;
// myfile.open("pos.csv");

// ofstream myfile2;
// myfile2.open("orientations.csv");

// ofstream myfile3;
// myfile3.open("patchpositions.csv");

// ofstream myfile4;
// myfile4.open("torques.csv");

// ofstream myfile5;
// myfile5.open("temperatures.csv");

int every = 100;


for(int i = 0 ; i < 500000 ; i++ ) {

//cout << i << endl;

a.advancemom_halfstep(F,T);
a.advance_pos();
a.rotate();

F = a.calculateforces(*pairs, wsa);
T = zeromatrix;
a.calculate_forces_and_torques3D(*pairs, pots, F, T);





a.create_forces_and_torques_sphere(F, T);



a.advancemom_halfstep(F,T);


// cout << a.getmom() << endl;
// cout << a.getangmom() << endl;
// pausel();
if(i%every== 0) {
    delete pairs;
    pairs = a.calculatepairs(boxes, 3.5);

    cout << i << endl;

    stringstream ss;
    ss << setw(5) << setfill('0') << (i / every);

    matrix<double> orient = a.getorientation();
    matrix<double> pos = a.getdat();

    string poss = "pos_i=";
    string oris = "orientation_i=";
    string extension =  ".csv";


    poss += ss.str();
    oris += ss.str();

    poss += extension;
    oris += extension; 

    ofstream myfile;
    myfile.open(poss.c_str());

    ofstream myfile2;
    myfile2.open(oris.c_str());
    // double qtemp0 = orient.gpcons(0, 0);
    // double qtemp1 = orient.gpcons(0, 1);
    // double qtemp2 = orient.gpcons(0, 2);
    // double qtemp3 = orient.gpcons(0, 3);
    // double qtemp4 = orient.gpcons(0, 4);
    // double qtemp5 = orient.gpcons(0, 5);
    // double qtemp6 = orient.gpcons(0, 6);
    // double qtemp7 = orient.gpcons(0, 7);
    // double qtemp8 = orient.gpcons(0, 8);

    // double gtemp0 = orient.gpcons(1, 0);
    // double gtemp1 = orient.gpcons(1, 1);
    // double gtemp2 = orient.gpcons(1, 2);
    // double gtemp3 = orient.gpcons(1, 3);
    // double gtemp4 = orient.gpcons(1, 4);
    // double gtemp5 = orient.gpcons(1, 5);
    // double gtemp6 = orient.gpcons(1, 6);
    // double gtemp7 = orient.gpcons(1, 7);
    // double gtemp8 = orient.gpcons(1, 8);

    // double nx1 = pot.nxb1 * qtemp0 + pot.nyb1 * qtemp3 + pot.nzb1 * qtemp6;
    // double ny1 = pot.nxb1 * qtemp1 + pot.nyb1 * qtemp4 + pot.nzb1 * qtemp7;
    // double nz1 = pot.nxb1 * qtemp2 + pot.nyb1 * qtemp5 + pot.nzb1 * qtemp8;

    // double nx2 = pot.nxb2 * gtemp0 + pot.nyb2 * gtemp3 + pot.nzb2 * gtemp6;
    // double ny2 = pot.nxb2 * gtemp1 + pot.nyb2 * gtemp4 + pot.nzb2 * gtemp7;
    // double nz2 = pot.nxb2 * gtemp2 + pot.nyb2 * gtemp5 + pot.nzb2 * gtemp8;

    // cout << pos(0, 0) + 0.5 * nx1 << "," << pos(0, 1) + 0.5 * ny1 << "," << pos(0, 2) + 0.5 * nz1 << endl;
    // cout << pos(1, 0) + 0.5 * nx2 << "," << pos(1, 1) + 0.5 * ny2 << "," << pos(1, 2) + 0.5 * nz2 << endl;
    // myfile3 << pos(0, 0) + 0.5 * nx1 << "," << pos(0, 1) + 0.5 * ny1 << "," << pos(0, 2) + 0.5 * nz1 << endl;
    // myfile3 << pos(1, 0) + 0.5 * nx2 << "," << pos(1, 1) + 0.5 * ny2 << "," << pos(1, 2) + 0.5 * nz2 << endl;

    myfile <<= pos;
    myfile2 <<= orient;

    myfile.close();
    myfile2.close();
    //myfile4 <<= tempT;

}

}
// myfile.close();
// myfile2.close();
// myfile3.close();
// myfile4.close();

*/

return 0;
}