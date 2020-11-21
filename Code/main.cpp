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

int n = 2000;

Condensate A(cbrt(2)*27.144, n);

//for one hundeed twenty five
// double basex =  10.0;
// double basey =  10.0;
// double basez =  10.0;

// double dx =1.0;

// int iter  = 0;
// matrix<double> initialpos(n,3);

// for(int i = 0  ; i < 5 ; i++) {
//     for(int j =  0 ; j < 5 ; j ++) {
//         for(int k = 0 ; k < 5 ; k++ ) {
//             initialpos(iter, 0) = basex + i * dx;
//             initialpos(iter, 1) = basey + j * dx;
//             initialpos(iter, 2) = basez + k * dx;
//             iter++;
//         }
//     }
// }

// A.obj->setdat(initialpos);

//TetrahedralPatch c(10.0, 1.4, pi / 4.);

TwoTetrahedral c(10.0, 1.4, pi / 4., 0.0, 1., pi / 6., 0.0, 1., pi / 6., 1000, 1000);

// int i = 2;

// int **p = new int*;
// int **q = new int*;

// *q = &i;

// *p = *q;

// cout << **q << endl;
// cout << **p << endl;

// int j =3 ;
// *q = &j;

// cout << **q << endl;
// cout << **p << endl;

// delete p;
// delete q;

// pausel();

// #pragma omp parallel for
// for(int i = 0  ; i < 10000 ; i++) {
// int p1 = rand() % 2000;
// int p2 = rand() % 2000;

// int **q = new int *;
// c.UpdateIteratorSafe(p1,p2,q);
// //*q = *(c.p);
// //c.assign_pointer(p1,p2,q);

// //cout << p1 << " " << p2 << " " << (*(c.p))[1] << endl
//     for (int tp = 1; tp < (*q)[0] + 1; tp++)
//     {

//     int potn = (*q)[tp];
//     // std::stringstream stream;
//     // stream << p1 << " " << p2 << " " << potn << endl;
//     if(p1 < 1000 && p2 < 1000 && potn > 15) {
//         //cout << stream.str();
//         error("found error");
//     }
//     else if(p1>1000 && p2 > 1000 && potn< 32) {
//         //cout << stream.str();
//         std::stringstream stream;
//         stream << p1 << " " << p2 << " " << potn << endl;
//         cout << stream.str();
//         error("found error");
        
//     }
//     else if(p1 > 1000 && p2 <1000 && (potn <16 || potn>31) ) {
//         //cout << stream.str();
//         error("found error");
//     }
//     else if (p1 < 1000 && p2 > 1000 && (potn < 16 || potn > 31))
//     {
//         //cout << stream.str();
//         error("found error");
//     }
//     else{

//     }
//     //cout << stream.str();
//     }
//     delete q;
// }

// //pausel();
// pausel(); 

BindingModelSingle b(0.99,0.01);

A.setBindingModel(b);

A.setpots(c);

int a = system("python3 /home/dino/Documents/Condensate/Code/Plotting/FigureMonitor.py ./ >filecreationlog &");

for (double kT = 1.0; kT > 0.49; kT -= 0.1)
{
    A.obj->setkT(kT);
    stringstream ss;
    ss << kT;

    string base = "_kT=";
    base += ss.str();

    A.run_singlebond(1000000, 1000, base);
}

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