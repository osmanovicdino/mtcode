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
int *p_indices_sorted = new int [N];

for(int i = 0 ; i < N ; i++)
p_indices[i]=i;
// float2 *d_particles = new float2 [N];
// float2 *d_cells = new float2 [N];


matrix<double> store(N,2);

for(int i = 0  ; i < N ; i++) {

float2 c;
c.x=5.*(rand()/(double)RAND_MAX);
c.y=5.*(rand()/(double)RAND_MAX);

store(i,0)=c.x;
store(i,1)=c.y;

(particles)[i]=c;
}


vector1<bool> pb(2,true);
cube bc(5.,pb,2);
int num = floor(5.0/1.0);
LangevinNVT test(bc);


test.setdat(store);
int ccc;
matrix<int> boxes = (test).getgeo().generate_boxes_relationships(4,ccc);
matrix<int> *froyo1 = test.calculatepairs(boxes,sqrt(1.75));

// cout << "CPU pairs" << endl;
// cout << *froyo1 << endl;

WCAPotential faa(1.0,1.0,1.0); 
matrix<double> F1(test.calculateforces(*froyo1,faa));
cout << "CPU forces" << endl;
cout << F1 << endl;



// cout << "particles: " << endl;
// for(int i =  0 ; i < N ; i++)
// cout << particles[i].x << " " << particles[i].y <<endl;

float2 *d_particles;
int *d_cells;
int *d_p_indices;





int size =  N*sizeof(float2);
int size2 = N*sizeof(int);


cudaMalloc((void**)&d_particles,size);
cudaMalloc((void**)&d_cells,size2);
cudaMalloc((void**)&d_p_indices,size2);

cudaMemcpy(d_particles,particles,size,cudaMemcpyHostToDevice);
cudaMemcpy(d_p_indices,p_indices,size2,cudaMemcpyHostToDevice);

assign_cell<<<N,1>>>(d_cells,d_particles,N,5.0,4.);

thrust::device_ptr<int> t_cells(d_cells);
thrust::device_ptr<int> t_indices(d_p_indices);

thrust::sort_by_key(t_cells,t_cells+N,t_indices);

thrust::device_vector<int>::iterator iter = thrust::max_element(t_cells,t_cells+N);

int largest = *iter;
// int largest = t_cells[iter];
// cout << "largest: " << largest << endl;

cudaMemcpy(cells,d_cells,size2,cudaMemcpyDeviceToHost);

cudaMemcpy(p_indices_sorted,d_p_indices,size2,cudaMemcpyDeviceToHost);

// cout << "cells: " << endl;
// for(int i =  0 ; i < N ; i++)
// cout << cells[i] << " " << p_indices_sorted[i] << endl;


int nperl = 4;
int ncells = nperl*nperl;
int *cellsc = new int [ncells];
for(int i = 0 ; i < ncells ; i++) {
	cellsc[i]=0;
}

int *d_cellsc;

int size3 = ncells*sizeof(int);
cudaMalloc((void**)&d_cellsc,size3);
cudaMemcpy(d_cellsc,cellsc,size3,cudaMemcpyHostToDevice);


cell_counts<<<N,1>>>(d_cells,d_cellsc,N,largest);


int *cellc2 = new int [ncells];

cudaMemcpy(cellc2,d_cellsc,size3,cudaMemcpyDeviceToHost);

// cout << "cell counts" << endl;
// for(int i =  0 ; i < ncells ; i++)
// cout << i << " " << cellc2[i] << endl;


//we now have the count of each cell list
int nbpairs = 5*ncells;

int *cells1 = new int [nbpairs];
int *cells2 = new int [nbpairs];
int *npb = new int [nbpairs];

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

		//+(0,0)
		//+(1,0)
		//+(-1,1)
		//+(0,1)
		//+(1,1)

	}
}

for(int i = 0 ; i < nbpairs ; i++) {
	npb[i] = 0;
}


// cout << "box checks" << endl;
// for(int i = 0 ; i < nbpairs ; i++) {
// cout << cells1[i] <<  " " << cells2[i] << endl;
// }

int size4 = nbpairs*sizeof(int);

int *d_cells1;
int *d_cells2;
int *d_npb;

cudaMalloc((void**)&d_cells1,size4);

cudaMalloc((void**)&d_cells2,size4);

cudaMalloc((void**)&d_npb,size4);




cudaMemcpy(d_cells1,cells1,size4,cudaMemcpyHostToDevice);
cudaMemcpy(d_cells2,cells2,size4,cudaMemcpyHostToDevice);
cudaMemcpy(d_npb,npb,size4,cudaMemcpyHostToDevice);

neighborlist_number<<<nbpairs,1>>>(d_cells1, d_cells2, d_cellsc,d_npb,nbpairs);


thrust::device_ptr<int> t_npb(d_npb);
//int tp =  thrust::reduce(t_npb,t_npb+nbpairs);

thrust::inclusive_scan(t_npb,t_npb+nbpairs,t_npb);

cudaMemcpy(npb,d_npb,size4,cudaMemcpyDeviceToHost);

// cout << "no.per box" << endl;
// for(int i = 0  ; i < nbpairs ; i++)
// 	cout << cells1[i] << " " <<cells2[i] << " " << npb[i] << endl;



// cout << "total" << " " << npb[nbpairs-1] << endl;

//int tpp =  npb[nbpairs-1];


int tpp;
cudaMemcpy(&tpp,d_npb+nbpairs-1,sizeof(int),cudaMemcpyDeviceToHost);
//cout << tpp << endl;
//cout << hTargetVariable << endl;



int *indices1 = new int [tpp];
int *indices2 = new int [tpp];
int *close =  new int [tpp];
for(int i = 0 ; i < tpp ; i++) {
	indices1[i]=0;
	indices2[i]=0;
	close[i]=0;
}

int *d_indices1;
int *d_indices2;
int *d_close;

int size5 = tpp*sizeof(int);

cudaMalloc((void**)&d_indices1,size5);

cudaMalloc((void**)&d_indices2,size5);

cudaMalloc((void**)&d_close,size5);

cudaMemcpy(d_indices1,indices1,size5,cudaMemcpyHostToDevice);
cudaMemcpy(d_indices2,indices2,size5,cudaMemcpyHostToDevice);
cudaMemcpy(d_close,close,size5,cudaMemcpyHostToDevice);


possible_neighborlist<<<nbpairs,1>>>(d_cells1, d_cells2, d_cellsc, d_p_indices, d_npb, nbpairs, d_indices1, d_indices2, d_close, d_particles,5., true,1.75);

cudaMemcpy(indices1,d_indices1,size5,cudaMemcpyDeviceToHost);
cudaMemcpy(indices2,d_indices2,size5,cudaMemcpyDeviceToHost);
cudaMemcpy(close,d_close,size5,cudaMemcpyDeviceToHost);

// for(int i = 0  ; i < tpp ; i++) {
// 	cout << indices1[i] <<  " " << indices2[i] << " " << close[i] << endl;
// }
// cout << "ok boomer" << endl;

thrust::device_ptr<int> t_close(d_close);
//int tp =  thrust::reduce(t_npb,t_npb+nbpairs);

thrust::inclusive_scan(t_close,t_close+tpp,t_close);

cudaMemcpy(close,d_close,size5,cudaMemcpyDeviceToHost);

// for(int i = 0  ; i < tpp ; i++) {
// 	cout << indices1[i] <<  " " << indices2[i] << " " << close[i] << endl;
// }

int th = close[tpp-1];
int *list1 =  new int [th];
int *list2 =  new int [th];
int *list3 =  new int [th];
int *list4 =  new int [th];

int *d_list1;
int *d_list2;
int *d_list3;
int *d_list4;

cudaMalloc((void**)&d_list1,th*sizeof(int));

cudaMalloc((void**)&d_list2,th*sizeof(int));

cudaMalloc((void**)&d_list3,th*sizeof(int));

cudaMalloc((void**)&d_list4,th*sizeof(int));

// cudaMemcpy(d_list1,list1,th*sizeof(int),cudaMemcpyHostToDevice);
// cudaMemcpy(d_list2,list2,th*sizeof(int),cudaMemcpyHostToDevice);

neighborlist<<<tpp,2>>>(d_indices1, d_indices2, d_close, tpp, d_list1, d_list2,d_list3,d_list4);

cudaMemcpy(list1,d_list1,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list2,d_list2,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list3,d_list3,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list4,d_list4,th*sizeof(int),cudaMemcpyDeviceToHost);



cout << "interacting species GPU" << endl;
for(int i = 0  ; i < th ; i++) {
	cout << list1[i] <<  " " << list2[i] << endl;
}

cout << "interacting species GPU" << endl;
for(int i = 0  ; i < th ; i++) {
	cout << list3[i] <<  " " << list4[i] << endl;
}

int yt;

int *d_p_indices2;
cudaMalloc((void**)&d_p_indices2,size2);
cudaMemcpy(d_p_indices2,p_indices,size2,cudaMemcpyHostToDevice);

int *d_list5;
int *d_list6;
int *d_list7;
int *d_list8;

construct_pair_list(d_particles,d_p_indices2,N,5.,d_cells1,d_cells2,4.,true,1.75,d_indices1,d_indices2,d_close);

cudaMemcpy(list1,d_list5,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list2,d_list6,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list3,d_list7,th*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(list4,d_list8,th*sizeof(int),cudaMemcpyDeviceToHost);



cout << "interacting species GPU from func" << endl;
for(int i = 0  ; i < th ; i++) {
	cout << list1[i] <<  " " << list2[i] << endl;
}

cout << "interacting species GPU from func" << endl;
for(int i = 0  ; i < th ; i++) {
	cout << list3[i] <<  " " << list4[i] << endl;
}


pausel();





double *d_forces1x;
double *d_forces1y;
double *d_forces2x;
double *d_forces2y;



cudaMalloc((void**)&d_forces1x,th*sizeof(double));


cudaMalloc((void**)&d_forces1y,th*sizeof(double));

cudaMalloc((void**)&d_forces2x,th*sizeof(double));


cudaMalloc((void**)&d_forces2y,th*sizeof(double));

gpupotential *inj = new gpupotential(1.0,1.0,1.0); 

gpupotential *d_inj;

cudaMalloc((void**)&d_inj,sizeof(gpupotential));

cudaMemcpy(d_inj,inj,sizeof(gpupotential),cudaMemcpyHostToDevice);


calculateforces2D<<<th,1>>>(d_list1,d_list2,d_particles,d_forces1x,d_forces1y,d_forces2x,d_forces2y, d_inj ,th, 5.,true);

double *forces1x = new double [th];

double *forces1y = new double [th];

double *forces2x = new double [th];

double *forces2y = new double [th];


// cudaMemcpy(forces1x,d_forces1x,th*sizeof(double),cudaMemcpyDeviceToHost);
// cudaMemcpy(forces2x,d_forces2x,th*sizeof(double),cudaMemcpyDeviceToHost);

// cudaMemcpy(forces1y,d_forces1y,th*sizeof(double),cudaMemcpyDeviceToHost);
// cudaMemcpy(forces2y,d_forces2y,th*sizeof(double),cudaMemcpyDeviceToHost);

// for(int i = 0 ; i < th ; i++) {
// 	cout << list1[i] << " " << list2[i] << " " <<forces1x[i] << " " << forces1y[i] << " " << forces2x[i] << " " << forces2y[i] << endl;
// }


thrust::device_ptr<double> t_forces1x(d_forces1x);
thrust::device_ptr<double> t_forces1y(d_forces1y);
thrust::device_ptr<int> t_list1(d_list1);
thrust::device_ptr<int> t_list3(d_list3);

thrust::device_ptr<double> t_forces2x(d_forces2x);
thrust::device_ptr<double> t_forces2y(d_forces2y);
thrust::device_ptr<int> t_list2(d_list2);
thrust::device_ptr<int> t_list4(d_list4);


thrust::sort_by_key(t_list1,t_list1+th,t_forces1x);


thrust::sort_by_key(t_list2,t_list2+th,t_forces2x);

thrust::sort_by_key(t_list3,t_list3+th,t_forces1y);

thrust::sort_by_key(t_list4,t_list4+th,t_forces2y);

// cudaMemcpy(list1,d_list1,th*sizeof(int),cudaMemcpyDeviceToHost);
// cudaMemcpy(list2,d_list2,th*sizeof(int),cudaMemcpyDeviceToHost);

// cudaMemcpy(forces1x,d_forces1x,th*sizeof(double),cudaMemcpyDeviceToHost);
// cudaMemcpy(forces2x,d_forces2x,th*sizeof(double),cudaMemcpyDeviceToHost);

// cudaMemcpy(forces1y,d_forces1y,th*sizeof(double),cudaMemcpyDeviceToHost);
// cudaMemcpy(forces2y,d_forces2y,th*sizeof(double),cudaMemcpyDeviceToHost);

// cout << "sorted" << endl;

// for(int i = 0 ; i < th ; i++) {
// 	cout << list1[i] << " " << list2[i] << " " <<forces1x[i] << " " << forces1y[i] << " " << forces2x[i] << " " << forces2y[i] << endl;
// }

double *d_sumforces1x;
double *d_sumforces1y;
double *d_sumforces2x;
double *d_sumforces2y;


int *d_key_reduce1x;
int *d_key_reduce1y;
int *d_key_reduce2x;
int *d_key_reduce2y;

cudaMalloc((void**)&d_sumforces1x,th*sizeof(double));

cudaMalloc((void**)&d_sumforces2x,th*sizeof(double));

cudaMalloc((void**)&d_sumforces1y,th*sizeof(double));

cudaMalloc((void**)&d_sumforces2y,th*sizeof(double));



cudaMalloc((void**)&d_key_reduce1x,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce2x,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce1y,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce2y,th*sizeof(int));


thrust::device_ptr<double> t_sumforces1x(d_sumforces1x);
thrust::device_ptr<double> t_sumforces2x(d_sumforces2x);
thrust::device_ptr<double> t_sumforces1y(d_sumforces1y);
thrust::device_ptr<double> t_sumforces2y(d_sumforces2y);

thrust::device_ptr<int> t_key_reduce1x(d_key_reduce1x);
thrust::device_ptr<int> t_key_reduce1y(d_key_reduce1y);
thrust::device_ptr<int> t_key_reduce2x(d_key_reduce2x);
thrust::device_ptr<int> t_key_reduce2y(d_key_reduce2y);


thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end1;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end2;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end3;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end4;

// cudaMemcpy(list1,d_list1,th*sizeof(int),cudaMemcpyDeviceToHost);
// cudaMemcpy(list2,d_list2,th*sizeof(int),cudaMemcpyDeviceToHost);


// cudaMemcpy(forces1x,d_forces1x,th*sizeof(double),cudaMemcpyDeviceToHost);
// cudaMemcpy(forces2x,d_forces2x,th*sizeof(double),cudaMemcpyDeviceToHost);

// cudaMemcpy(forces1y,d_forces1y,th*sizeof(double),cudaMemcpyDeviceToHost);
// cudaMemcpy(forces2y,d_forces2y,th*sizeof(double),cudaMemcpyDeviceToHost);

// for(int i = 0 ; i < th ; i++) {
// 	cout << list1[i] << " " << forces1x[i] << endl;
// }



new_end1 = thrust::reduce_by_key(t_list1,t_list1+th,t_forces1x,t_key_reduce1x,t_sumforces1x);


new_end2 = thrust::reduce_by_key(t_list3,t_list3+th,t_forces1y,t_key_reduce1y,t_sumforces1y);


new_end3 = thrust::reduce_by_key(t_list2,t_list2+th,t_forces2x,t_key_reduce2x,t_sumforces2x);


new_end4 = thrust::reduce_by_key(t_list4,t_list4+th,t_forces2y,t_key_reduce2y,t_sumforces2y);



// cout << "first: " << (new_end1.first) << endl;
// cout << "second: " << (new_end1.second) << endl;

// // int *raw_ptr = thrust::raw_pointer_cast(new_end1.first);
//cout << "hmm: " << thrust::raw_pointer_cast(&new_end1.second[0])-thrust::raw_pointer_cast(&t_sumforces1x[0]) << endl;
int ih1 = thrust::raw_pointer_cast(&new_end1.first[0])-thrust::raw_pointer_cast(&t_key_reduce1x[0]);
//int ih2 = thrust::raw_pointer_cast(&new_end2.first[0])-thrust::raw_pointer_cast(&t_key_reduce1y[0]);
int ih3 = thrust::raw_pointer_cast(&new_end3.first[0])-thrust::raw_pointer_cast(&t_key_reduce2x[0]);
//int ih4 = thrust::raw_pointer_cast(&new_end4.first[0])-thrust::raw_pointer_cast(&t_key_reduce2y[0]);
// cout << &raw_ptr[0] << endl;

//the output of this process will look like 

// cudaMemcpy(forces1x,d_sumforces1x,th*sizeof(double),cudaMemcpyDeviceToHost);


// for(int i = 0 ; i < th ; i++) {
// 	cout << list1[i] << " "  << forces1x[i] << endl;
// }

double *totalforcex = new double [N];
double *totalforcey = new double [N];

for(int i = 0  ; i < N ; i++) {
	totalforcex[i]=0.0;
	totalforcey[i]=0.0;
}


// cout << "before forces" << endl;
// for(int i = 0 ; i < N ; i++) {
// 	cout << totalforcex[i] << " " << totalforcey[i] << endl;
// }

double *d_totalforcex;
double *d_totalforcey;

cudaMalloc((void**)&d_totalforcex,N*sizeof(double));

cudaMalloc((void**)&d_totalforcey,N*sizeof(double));

cudaMemcpy(d_totalforcex,totalforcex,N*sizeof(double),cudaMemcpyHostToDevice);

cudaMemcpy(d_totalforcey,totalforcey,N*sizeof(double),cudaMemcpyHostToDevice);


addforce<<<ih1,1>>>(d_totalforcex,d_totalforcey, d_key_reduce1x,d_key_reduce1y, d_sumforces1x, d_sumforces1y, ih1);

addforce<<<ih3,1>>>(d_totalforcex,d_totalforcey, d_key_reduce2x,d_key_reduce2y, d_sumforces2x, d_sumforces2y, ih3);


cudaMemcpy(totalforcex,d_totalforcex,N*sizeof(double),cudaMemcpyDeviceToHost);

cudaMemcpy(totalforcey,d_totalforcey,N*sizeof(double),cudaMemcpyDeviceToHost);


cout << "reduced" << endl;

cout << "GPU forces" << endl;
for(int i = 0 ; i < N ; i++) {
	cout << totalforcex[i] << " " << totalforcey[i] << endl;
}


free(particles); free(cells);

cudaFree(d_particles); cudaFree(d_cells);







return 0;
}