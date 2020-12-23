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
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
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

// #include "NCGasR.h"
// #include "Microtubule.h"



//#include "MDGPU.cu"

using namespace std;



int main(int argc, char** argv) {

srand (time(NULL));

/*
Condensate A(4.0,2);

matrix<double> dat(2,3);

dat(0, 0) = 2.0;
dat(0, 1) = 2.0;
dat(0, 2) = 2.0;

dat(1,0) = 2.0;
dat(1,1) = 2.0;
dat(1,2) = 3.4;

A.obj->setdat(dat);

TetrahedralPatch c(10.0, 1.0, pi / 3.);
A.setpots(c);

vector1<double> v1(3);

v1[0] = 0.;
v1[1] = 0.;
v1[2] = 1.;

vector1<double> v12(3);

double dy = 0.6;
double dz = 0.6;

v12[0] = dy;
v12[1] = dz;
v12[2] = sqrt(SQR(1.) - SQR(dy) - SQR(dz));

vector1<double> v2(3);

v2[0] = -sqrt(2. / 9.);
v2[1] = sqrt(2. / 3.);
v2[2] = -1. / 3.;

vector1<double> v3(3);

v3[0] = -dy;
v3[1] = -dz;
v3[2] = -sqrt(SQR(1.) - SQR(dy) - SQR(dz));

A.obj->setorientation(v1, v12, 0);
A.obj->setorientation(v2, v3, 1);

// cout << A.obj->getorientation() << endl;

BindingModelSingle b(0.99, 0.00);

A.setBindingModel(b);

A.setviscosity(0.0);
A.setkT(0.0);

A.run(1000,1); 


*/
/* 
BindingModelBinary b(4*125);




//particles type 1 and 2 can bind to each other;

vector1<double> r11(4);
r11[0] = 0.0;
r11[1] = 1.0;
r11[2] = 0.0;
r11[3] = 1.0;

vector1<double> r12(r11);

vector1<double> r22(4);
r22[0] = 1.0;
r22[1] = 0.0;
r22[2] = 1.0;
r22[3] = 0.0;

vector1<double> t111(16);
vector1<double> t112(16);
vector1<double> t122(16);
vector1<double> t222(16);



t111[0]=0.998;//fromi,jtpi,j
t111[1]=0.001;//fromi,jtoj,k
t111[2]=0.001;//fromi,jtoi,k
t111[3]=0.000;//fromi,jtonothing
t111[4]=0.001;//fromj,ktoi,j
t111[5]=0.998;//fromj,ktoj,k
t111[6]=0.001;//fromj,ktoi,k
t111[7]=0.000;//fromj,ktonothing
t111[8]=0.001;//fromi,ktoi,j
t111[9]=0.001;//fromi,ktoj,k
t111[10]=0.998;//fromi,ktoi,k
t111[11]=0.000;//fromi,ktonothing
t111[12]=0.33;//fromnothingtoi,j
t111[13]=0.33;//fromnothingtoj,k
t111[14]=0.33;//fromnothingtoi,k
t111[15]=0.01;//fromnothingtonothing

//particle 2 does not interact

t112[0] = 0.001;  //fromi,jtpi,j
t112[1] = 0.499;  //fromi,jtoj,k
t112[2] = 0.499;  //fromi,jtoi,k
t112[3] = 0.000;  //fromi,jtonothing
t112[4] = 0.000;  //fromj,ktoi,j
t112[5] = 0.998;  //fromj,ktoj,k
t112[6] = 0.100;  //fromj,ktoi,k
t112[7] = 0.001;  //fromj,ktonothing
t112[8] = 0.000;  //fromi,ktoi,j
t112[9] = 0.100;  //fromi,ktoj,k
t112[10] = 0.998; //fromi,ktoi,k
t112[11] = 0.001; //fromi,ktonothing
t112[12] = 0.001; //fromnothingtoi,j
t112[13] = 0.499; //fromnothingtoj,k
t112[14] = 0.499; //fromnothingtoi,k
t112[15] = 0.000; //fromnothingtonothing

t122[0] = 0.998;  //fromi,jtpi,j
t122[1] = 0.000;  //fromi,jtoj,k
t122[2] = 0.100;  //fromi,jtoi,k
t122[3] = 0.001;  //fromi,jtonothing
t122[4] = 0.499;  //fromj,ktoi,j
t122[5] = 0.000;  //fromj,ktoj,k
t122[6] = 0.499;  //fromj,ktoi,k
t122[7] = 0.001;  //fromj,ktonothing
t122[8] = 0.100;  //fromi,ktoi,j
t122[9] = 0.000;  //fromi,ktoj,k
t122[10] = 0.998; //fromi,ktoi,k
t122[11] = 0.001; //fromi,ktonothing
t122[12] = 0.499; //fromnothingtoi,j
t122[13] = 0.000; //fromnothingtoj,k
t122[14] = 0.499; //fromnothingtoi,k
t122[15] = 0.001; //fromnothingtonothing

t222[0] = 0.000;  //fromi,jtpi,j
t222[1] = 0.000;  //fromi,jtoj,k
t222[2] = 0.000;  //fromi,jtoi,k
t222[3] = 1.0;    //fromi,jtonothing
t222[4] = 0.000;  //fromj,ktoi,j
t222[5] = 0.000;  //fromj,ktoj,k
t222[6] = 0.000;  //fromj,ktoi,k
t222[7] = 1.000;  //fromj,ktonothing
t222[8] = 0.000;  //fromi,ktoi,j
t222[9] = 0.000;  //fromi,ktoj,k
t222[10] = 0.000; //fromi,ktoi,k
t222[11] = 1.000; //fromi,ktonothing
t222[12] = 0.000; //fromnothingtoi,j
t222[13] = 0.000; //fromnothingtoj,k
t222[14] = 0.000; //fromnothingtoi,k
t222[15] = 1.000; //fromnothingtonothing

b.doubr11 = r11;
b.doubr12 = r12;
b.doubr22 = r22;

b.tripr111 = t111;
b.tripr112 = t112;
b.tripr122 = t122;
b.tripr222 = t222;

TetrahedralWithSingle c(30.0,1.4,pi/4.,100.,1.4,pi/3.,0.0,1.,pi/6.,125 , 400);

double T;
bool err1,err2;
matrix<double> initpos = importcsv("./Basic/InitialConditions/initpos.csv", T, err1);
matrix<double> initori = importcsv("./Basic/InitialConditions/initori.csv", T, err2);





Condensate A(21.5443,525);

A.obj->setdat(initpos);
A.obj->setorientation(initori);

A.setBindingModel(b);

A.setpots(c);



A.run_singlebond(1000000, 1000);
 */

/*
vector1<int> indexes(1000);
vector1<int> indexes_np(1000);
int T;
bool er;

matrix<int> pairs = importcsv("./Basic/InitialConditions/pairs.csv",T,er);

matrix<int> possibles(1000, 100);
matrix<int> possibles_np(1000, 100);


std::mutex mtx;

#pragma omp parallel for
for(int i = 0 ; i < pairs.getNsafe() ; i++ ) {
    int i1 =  pairs(i,0);
    int i2 =  pairs(i,1);


    // #pragma atomic read 
    // int iterator1 = indexes[i1];
    // #pragma atomic read
    // int iterator2 = indexes[i2];
  

    mtx.lock();
    int iterator1 = indexes[i1];
    int iterator2 = indexes[i2];

    indexes[i1]++;
    indexes[i2]++;
    mtx.unlock();


    possibles(i1, iterator1) = i2;
    possibles(i2, iterator2) = i1;

    // #pragma atomic write
    // possibles(i2, temp.b) = i1;



    //    indexes[i2]++;
    

}



for (int i = 0; i < pairs.getNsafe(); i++)
{
    int i1 = pairs(i, 0);
    int i2 = pairs(i, 1);


    possibles_np(i1, indexes_np[i1]) = i2;
    possibles_np(i2, indexes_np[i2]) = i1;

    indexes_np[i1]++;
    indexes_np[i2]++;
}



cout << (indexes == indexes_np) << endl;

// cout << possibles(0,'r') << endl;
// cout << possibles_np(0,'r') << endl;



for(int i = 0  ; i < possibles.getNsafe() ; i++) {
    if(meanish(possibles_np(i,'r')) == meanish(possibles(i,'r'))) {

    }
    else {
        cout << possibles(i, 'r') << endl;
        cout << possibles_np(i, 'r') << endl;
        pausel();
    }
}
cout << "done" << endl;
pausel();
*/



BindingModelTernary b(2000*4,4*4000);


// b.setup(0.99,0.01,0.01,0.01,0.0,0.0,
// 0.,
// 0.0,
// 0.0,
// 0.0,
// 0.0,
// 0.0,
// 0.0,
// 0.0,
// 0.0);


b.setup(0.99,0.01,0.01,0.01,0.99,0.2,
0.,
8.,
0.0,
0.0,
0.0,
1.0,
0.0,
0.0,
0.0);


matrix<double> params(6,3);
params(0, 0) = 10.0; // 1<->1
params(0, 1) = 1.4;
params(0, 2) = 0.927;

params(1, 0) = 0.0; //1 <->2
params(1, 1) = 1.0;
params(1, 2) = 0.927;

params(2, 0) = 10.0; //1 <-> 3
params(2, 1) = 1.4;
params(2, 2) = 0.927;

params(3, 0) = 0.0; // 2 <-> 2
params(3, 1) = 1.0;
params(3, 2) = 0.927;

params(4, 0) = 10.0; // 2<-> 3
params(4, 1) = 1.4;
params(4, 2) = 0.927;

params(5, 0) = 0.0; //3 <-> 3
params(5, 1) = 1.0;
params(5, 2) = 0.927;

TwoTetrahedralAndSingle c(params, 2000, 4000, 8000);





int n = 8000;
int n2 = 100;
double packing_fraction = 0.05;

double l = cbrt(pi*(double)n/(6.*packing_fraction));

Condensate A(l, n);


//TwoTetrahedral c(10.0, 1.4, pi / 4., 0.0, 1., pi / 6., 0.0, 1., pi / 6., 1000, 1000);

// string filp = "/home/dino/Desktop/Chemistry/SimulationResults/ChemicalOscillator/sim-20-12-14-19:43:58/pos_beta=1_i=0455.csv";
// string filo = "/home/dino/Desktop/Chemistry/SimulationResults/ChemicalOscillator/sim-20-12-14-19:43:58/orientation_beta=1_i=0455.csv";

// string filp = "/home/dino/Documents/Condensate/TernaryFluid2/pos_beta=1_i=02097.csv";
// string filo = "/home/dino/Documents/Condensate/TernaryFluid2/orientation_beta=1_i=02097.csv";

// double T;
// bool err1;
// bool err2;
// matrix<double> temppos = importcsv(filp, T, err1);
// matrix<double> tempori = importcsv(filo, T, err1);

// matrix<double> newpos(2000,3);
// matrix<double> newori(2000,3);

// A.obj->setdat(temppos);
// A.obj->setorientation(tempori);


// cout << b.tripr111 << endl;
// cout << b.doubr11 << endl;
// cout << b.doubr22 << endl;
// cout << b2.on_rate << endl;
// cout << b2.off_rate << endl;
// pausel();

//TetrahedralPatch c2(10.0,1.4,0.927);

A.setBindingModel(b);

A.setpots(c);

//int a = system("python3 /home/dino/Documents/Condensate/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

//int a = system("python3 /home/dino/Desktop/tylercollab/Repo/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

//A.run_singlebond(10000, 1000);

A.setviscosity(1.0);


double beta = 1.0;

A.obj->setkT(1./beta);



stringstream ss;
ss << beta;

string base = "_beta=";
base += ss.str();

A.run_singlebond(10000000, 1000, base);


/* 
int ccc;
matrix<int> boxes = A.obj->getgeo().generate_boxes_relationships(A.num, ccc);
matrix<int> *pairs = A.obj->calculatepairs(boxes, 3.5);

int NN = A.obj->getN();
BinaryBindStore bbs;

int nh = (A.pots)->get_total_patches(NN);

vector1<bool> isbound(nh);

vector1<int> boundto(nh);

bbs.isbound = isbound;
bbs.boundto = boundto;

matrix<double> F(NN, 3);
matrix<double> F2(NN, 3);
matrix<double> T(NN, 3);


A.obj->calculate_forces_and_torques3D_onlyone(*pairs, c2, bbs, b2, F,T);

//cout << F << endl;
//cout << T << endl;

vector1<int> ty =  bbs.boundto;
//F.reset(0.0);
T.reset(0.0);

bbs.isbound = isbound;
bbs.boundto = boundto;

A.obj->calculate_forces_and_torques3D_onlyone(*pairs, c, bbs, b, F2, T);

// pausel();

// for(int i = 0  ; i < NN ; i++) {
//     cout << ty[i] << " " << bbs.boundto[i] << endl;
// }



for(int i = 0  ; i < NN ; i++) {
    cout << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << "\t" << F2(i, 0) << " " << F2(i, 1) << " " << F2(i, 2) << endl;
}
 */

// pausel();

//cout << T << endl;
//A.run_singlebond(1000000, 1000, base);

/* 
int n = 4;

Condensate A(5.0, n); //Create a condensate


BindingModelFull b(n);


matrix<double> matb(n*n,4);

for(int i = 0  ; i < n ; i++) {
    for(int j = i+1  ; j < n ; j++) {
        int index1 = i*n+j;
        if(i == 0) { //we define particle 0 as the promiscuous particle (can bind to anything)
            matb(index1,0) = 0.0; //unbound->unbound
            matb(index1,1) = 1.0; //unbound->bound
            matb(index1,2) = 0.0; //bound -> unbound
            matb(index1,3) = 1.0; //bound -> bound
        }
        else{ //no other particles can bind
            matb(index1, 0) = 1.0; //unbound->unbound
            matb(index1, 1) = 0.0; //unbound->bound
            matb(index1, 2) = 1.0; //bound -> unbound
            matb(index1, 3) = 0.0; //bound -> bound
        }
    }
}





b.doubrates = matb;



matrix<double> mat1(4, 4);
matrix<double> mat2(4, 4);
matrix<double> mat3(4, 4);
matrix<double> mat4(4, 4);


//i ==0  j ==1 k == 2
mat1(0,0) = 0.001;
mat1(0,1) = 0.001;
mat1(0,2) = 0.998;
mat1(0,3) = 0.0;
mat1(1,0) = 0.001;
mat1(1,1) = 0.001;
mat1(1,2) = 0.098;
mat1(1,3) = 0.00;
mat1(2,0) = 0.001;
mat1(2,1) = 0.001;
mat1(2,2) = 0.998;
mat1(2,3) = 0.0;
mat1(3,0) = 0.001;
mat1(3,1) = 0.001;
mat1(3,2) = 0.998;
mat1(3,3) = 0.0;

//i ==0  j ==1 k == 3
mat2(0, 0) = 0.998;
mat2(0, 1) = 0.001;
mat2(0, 2) = 0.001;
mat2(0, 3) = 0.0;
mat2(1, 0) = 0.998;
mat2(1, 1) = 0.001;
mat2(1, 2) = 0.001;
mat2(1, 3) = 0.00;
mat2(2, 0) = 0.998;
mat2(2, 1) = 0.001;
mat2(2, 2) = 0.001;
mat2(2, 3) = 0.0;
mat2(3, 0) = 0.998;
mat2(3, 1) = 0.001;
mat2(3, 2) = 0.001;
mat2(3, 3) = 0.0;


//i ==0  j ==2 k == 3
mat3(0, 0) = 0.001;
mat3(0, 1) = 0.001;
mat3(0, 2) = 0.998;
mat3(0, 3) = 0.0;
mat3(1, 0) = 0.001;
mat3(1, 1) = 0.001;
mat3(1, 2) = 0.998;
mat3(1, 3) = 0.00;
mat3(2, 0) = 0.001;
mat3(2, 1) = 0.001;
mat3(2, 2) = 0.998;
mat3(2, 3) = 0.0;
mat3(3, 0) = 0.001;
mat3(3, 1) = 0.001;
mat3(3, 2) = 0.998;
mat3(3, 3) = 0.0;

//i ==1  j ==2 k == 3
mat4(0, 0) = 0.0;
mat4(0, 1) = 0.0;
mat4(0, 2) = 0.0;
mat4(0, 3) = 1.0;
mat4(1, 0) = 0.0;
mat4(1, 1) = 0.0;
mat4(1, 2) = 0.0;
mat4(1, 3) = 1.0;
mat4(2, 0) = 0.0;
mat4(2, 1) = 0.0;
mat4(2, 2) = 0.0;
mat4(2, 3) = 1.0;
mat4(3, 0) = 0.0;
mat4(3, 1) = 0.0;
mat4(3, 2) = 0.0;
mat4(3, 3) = 1.0;

matrix<double> triprates(n*n*n,16);
for(int i  = 0  ; i < n ; i++) {
    for(int j = i+1 ; j < n ; j++) {
        for(int k = j+1 ; k < n ; k++) {
          //  cout << i << " " << j << " " << k << endl;
            int indx = i*n*n+j*n+k;
            //cout << indx << endl;
            if(i == 0 && j ==1 && k == 2) {
                triprates(indx, 0) =  mat1(0, 0); //from i,j tp i,j
                triprates(indx, 1) =  mat1(0, 1);  //from i,j to j,k
                triprates(indx, 2) =  mat1(0, 2);  //from i,j to i,k
                triprates(indx, 3) =  mat1(0, 3);  //from i,j to nothing
                triprates(indx, 4) =  mat1(1, 0);  //from j,k to i,j
                triprates(indx, 5) =  mat1(1, 1);  //from j,k to j,k
                triprates(indx, 6) =  mat1(1, 2);  //from j,k to i,k
                triprates(indx, 7) =  mat1(1, 3);  //from j,k to nothing
                triprates(indx, 8) =  mat1(2, 0);  //from i,k to i,j
                triprates(indx, 9) =  mat1(2, 1);  //from i,k to j,k
                triprates(indx, 10) = mat1(2, 2); //from i,k to i,k
                triprates(indx, 11) = mat1(2, 3); //from i,k to nothing
                triprates(indx, 12) = mat1(3, 0); //from nothing to i,j
                triprates(indx, 13) = mat1(3, 1); //from nothing to j,k
                triprates(indx, 14) = mat1(3, 2); //from nothing to i,k
                triprates(indx, 15) = mat1(3, 3); //from nothing to nothing
            }
            else if(i==0 && j==1 && k ==3) 
            {
                triprates(indx, 0) = mat2(0, 0);  //from i,j tp i,j
                triprates(indx, 1) = mat2(0, 1);  //from i,j to j,k
                triprates(indx, 2) = mat2(0, 2);  //from i,j to i,k
                triprates(indx, 3) = mat2(0, 3);  //from i,j to nothing
                triprates(indx, 4) = mat2(1, 0);  //from j,k to i,j
                triprates(indx, 5) = mat2(1, 1);  //from j,k to j,k
                triprates(indx, 6) = mat2(1, 2);  //from j,k to i,k
                triprates(indx, 7) = mat2(1, 3);  //from j,k to nothing
                triprates(indx, 8) = mat2(2, 0);  //from i,k to i,j
                triprates(indx, 9) = mat2(2, 1);  //from i,k to j,k
                triprates(indx, 10) = mat2(2, 2); //from i,k to i,k
                triprates(indx, 11) = mat2(2, 3); //from i,k to nothing
                triprates(indx, 12) = mat2(3, 0); //from nothing to i,j
                triprates(indx, 13) = mat2(3, 1); //from nothing to j,k
                triprates(indx, 14) = mat2(3, 2); //from nothing to i,k
                triprates(indx, 15) = mat2(3, 3); //from nothing to nothing
            }
            else if (i == 0 && j == 2 && k == 3)
            {
                triprates(indx, 0) = mat3(0, 0);  //from i,j tp i,j
                triprates(indx, 1) = mat3(0, 1);  //from i,j to j,k
                triprates(indx, 2) = mat3(0, 2);  //from i,j to i,k
                triprates(indx, 3) = mat3(0, 3);  //from i,j to nothing
                triprates(indx, 4) = mat3(1, 0);  //from j,k to i,j
                triprates(indx, 5) = mat3(1, 1);  //from j,k to j,k
                triprates(indx, 6) = mat3(1, 2);  //from j,k to i,k
                triprates(indx, 7) = mat3(1, 3);  //from j,k to nothing
                triprates(indx, 8) = mat3(2, 0);  //from i,k to i,j
                triprates(indx, 9) = mat3(2, 1);  //from i,k to j,k
                triprates(indx, 10) = mat3(2, 2); //from i,k to i,k
                triprates(indx, 11) = mat3(2, 3); //from i,k to nothing
                triprates(indx, 12) = mat3(3, 0); //from nothing to i,j
                triprates(indx, 13) = mat3(3, 1); //from nothing to j,k
                triprates(indx, 14) = mat3(3, 2); //from nothing to i,k
                triprates(indx, 15) = mat3(3, 3); //from nothing to nothing
            }
            else if (i == 1 && j == 2 && k == 3)
            {
                triprates(indx, 0) = mat4(0, 0);  //from i,j tp i,j
                triprates(indx, 1) = mat4(0, 1);  //from i,j to j,k
                triprates(indx, 2) = mat4(0, 2);  //from i,j to i,k
                triprates(indx, 3) = mat4(0, 3);  //from i,j to nothing
                triprates(indx, 4) = mat4(1, 0);  //from j,k to i,j
                triprates(indx, 5) = mat4(1, 1);  //from j,k to j,k
                triprates(indx, 6) = mat4(1, 2);  //from j,k to i,k
                triprates(indx, 7) = mat4(1, 3);  //from j,k to nothing
                triprates(indx, 8) = mat4(2, 0);  //from i,k to i,j
                triprates(indx, 9) = mat4(2, 1);  //from i,k to j,k
                triprates(indx, 10) = mat4(2, 2); //from i,k to i,k
                triprates(indx, 11) = mat4(2, 3); //from i,k to nothing
                triprates(indx, 12) = mat4(3, 0); //from nothing to i,j
                triprates(indx, 13) = mat4(3, 1); //from nothing to j,k
                triprates(indx, 14) = mat4(3, 2); //from nothing to i,k
                triprates(indx, 15) = mat4(3, 3); //from nothing to nothing
            }
            else{
                //
            }
        }
    }
}



b.triprates =  triprates;

SingPatch c(100.0,2.,pi/3.);

A.setBindingModel(b);

A.setpots(c);

A.obj->setkT(1.0);


A.run_singlebond( 100000,  100); */
//SAVE CC
// vector1<bool> already_accounted(200, false);
// vector1<int> indexes(200,0);
// int iter = 0;
// vector1<int> demarkus(200,0);
// int len_dem = 0;

// int num_con_comp = 0;

// for (int i = 0; i < 200; i++) {

//     if (already_accounted[i] == false)
//     {
//         DFUtilstore(i, already_accounted, adj, lens, indexes, iter, len_dem);

//         demarkus[num_con_comp] = len_dem;

//         num_con_comp++;
//         len_dem = 0;

//        // cout << "\n";
//     }
// }

// cout << indexes << endl;
// cout << demarkus << endl;

//PRINT CC
// for (int i = 0; i < 200; i++)
// {

//     if (already_accounted[i] == false)
//     {
//         DFUtil(i, already_accounted, adj, lens);

//         // demarkus[num_con_comp] = len_dem;

//         // num_con_comp++;
//         // len_dem = 0;

//          cout << "\n";
//     }
// }
// Condensate A(25.0, 1000, 2);

// vector1<double> as(16,10.0);
// vector1<double> bs(16,1.5);
// vector1<double> cs(16,pi/3.);

// //as[0] = 100.0;

// set_potential_bundle_tetrahedral(A, as, bs, cs);

// for(double kT = 1.0 ; kT > 0.49 ; kT -= 0.1) {

// A.obj->setkT(kT);
// stringstream ss;
// ss << kT;

// string base = "_kT=";
// base += ss.str();

// A.run(1000000,1000,base);
// }

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