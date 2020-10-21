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
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "basic.h"
#include "vector1.h"
#include "matrix2.h"
#include "matrix2.cpp"
#include "potential.h"
#include "MD.h"
#include "Langevin.h"


// #include "BrownianGel.cpp"
// #include "BrownianGel2.cpp"
// #include "LangevinGel.cpp"
// #include "LangevinGelFixed.cpp"

// #include "NCGasR.h"
//#include "Microtubule.h"


#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include "MDGPU.cu"

using namespace std;



void prdshft(int &i, int max) {
if(i >= max) i=i-max;
else if(i<0) i=i+max;
else {}
}



int main(int argc, char** argv) {
srand (time(NULL));






int N = 10;
float2 *particles = new float2 [N];
int *cells = new int [N];
int *p_indices = new int [N];

for(int i = 0 ; i < N ; i++)
p_indices[i]=i;
// float2 *d_particles = new float2 [N];
// float2 *d_cells = new float2 [N];


for(int i = 0  ; i < N ; i++) {

float2 c;
c.x=5.*(rand()/(double)RAND_MAX);
c.y=5.*(rand()/(double)RAND_MAX);

(particles)[i]=c;
}




// cout << "particles: " << endl;
// for(int i =  0 ; i < N ; i++)
// cout << particles[i].x << " " << particles[i].y <<endl;

float2 *d_particles;
int *d_p_indices;





int size =  N*sizeof(float2);
int size2 = N*sizeof(int);


cudaMalloc((void**)&d_particles,size);
cudaMalloc((void**)&d_p_indices,size2);

cudaMemcpy(d_particles,particles,size,cudaMemcpyHostToDevice);
cudaMemcpy(d_p_indices,p_indices,size2,cudaMemcpyHostToDevice);

// cout << "cells: " << endl;
// for(int i =  0 ; i < N ; i++)
// cout << cells[i] << " " << p_indices_sorted[i] << endl;


int nperl = 4;
int ncells = nperl*nperl;


// cout << "cell counts" << endl;
// for(int i =  0 ; i < ncells ; i++)
// cout << i << " " << cellc2[i] << endl;


//we now have the count of each cell list
int nbpairs = 5*ncells;

int *cells1 = new int [nbpairs];
int *cells2 = new int [nbpairs];


int itery = 0;
for(int i1 = 0 ; i1 < nperl ; i1++) {
	for(int i2 = 0 ; i2 < nperl ; i2++ ) {


		int b1 =  i1*nperl+i2;

		int i3 = i1+0;
		int j3 = i2+0;

		int i4 = i1+1;
		int j4 = i2+0;

		int i5 = i1-1;
		int j5 = i2+1;

		int i6 = i1+0;
		int j6 = i2+1;

		int i7 = i1+1;
		int j7 = i2+1;

		prdshft(i3,nperl);
		prdshft(j3,nperl);

		prdshft(i4,nperl);
		prdshft(j4,nperl);

		prdshft(i5,nperl);
		prdshft(j5,nperl);
		
		prdshft(i6,nperl);
		prdshft(j6,nperl);
		
		prdshft(i7,nperl);
		prdshft(j7,nperl);		

		cells1[itery] =  b1;
		cells2[itery] =  i3*nperl+j3;

		itery++;

		cells1[itery] =  b1;
		cells2[itery] =  i4*nperl+j4;
		
		itery++;
		
		cells1[itery] =  b1;
		cells2[itery] =  i5*nperl+j5;
		
		itery++;
		
		cells1[itery] =  b1;
		cells2[itery] =  i6*nperl+j6;
		
		itery++;
		
		cells1[itery] =  b1;
		cells2[itery] =  i7*nperl+j7;

		itery++;


	}
}
int size4 = nbpairs*sizeof(int);

int *d_cells1;
int *d_cells2;

cudaMalloc((void**)&d_cells1,size4);

cudaMalloc((void**)&d_cells2,size4);

cudaMemcpy(d_cells1,cells1,size4,cudaMemcpyHostToDevice);
cudaMemcpy(d_cells2,cells2,size4,cudaMemcpyHostToDevice);

int *d_indices1;
int *d_indices2;
double *d_close;


int tpp;

construct_possible_pair_list(d_particles,d_p_indices,N,5.,d_cells1,d_cells2,4.,true,d_indices1,d_indices2,d_close,tpp);

// int *indices1 = new int [tpp];
// int *indices2 = new int [tpp];
// double *close = new double [tpp];
// cudaMemcpy(indices1,d_indices1,tpp*sizeof(int),cudaMemcpyDeviceToHost);
// cudaMemcpy(indices2,d_indices2,tpp*sizeof(int),cudaMemcpyDeviceToHost);
// cudaMemcpy(close,d_close,tpp*sizeof(double),cudaMemcpyDeviceToHost);

// for(int i = 0 ; i < tpp ; i++) {
// 	cout << indices1[i] << " " << indices2[i] << " " << close[i] << endl;
// }

// pausel();

//cout << tpp << endl;

int *d_list1;
int *d_list2;
int *d_list3;
int *d_list4;

less_than_condition_NAND cond1(1.75,0,5,5,10);


int th;
pairlist(d_indices1,d_indices2,d_close,cond1,d_list1,d_list2,d_list3,d_list4,tpp, th);





//cudaFree(d_list1);
cout << tpp << endl;
cout << th << endl;

int *list1 = new int [th];
int *list2 = new int [th];
int *list3 = new int [th];
int *list4 = new int [th];

cudaMemcpy(list1,d_list1,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list2,d_list2,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list3,d_list3,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list4,d_list4,th*sizeof(int),cudaMemcpyDeviceToHost);

for(int i = 0 ; i < th ; i++) {
	cout << list1[i] << " " << list2[i] << " " << list3[i] << " " << list4[i] << endl;
}

//pairlist(d_indices1,d_indices2,d_close,cond1,d_list1,d_list2,d_list3,d_list4);

// cudaMemcpy(p_indices,d_p_indices,size2,cudaMemcpyDeviceToHost);

// cout << "after gold" << endl;

// for(int i = 0 ; i < N ; i++)
// cout << p_indices[i] << endl;


return 0;
}