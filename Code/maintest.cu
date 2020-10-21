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
#include <thrust/unique.h>
#include <thrust/device_delete.h>
#include <curand.h>
#include <curand_kernel.h>
#include "MDGPU.cu"

using namespace std;

void create_list(int *&array) {

	int n = 5;
	cudaMalloc((void**)&array,n*sizeof(int));
	int h = 1;
	setstate<<<n,1>>>(array,1,5);

}

int main(int argc, char** argv) {

int *d_list;

create_list(d_list);

print_device_array(d_list,5);

cudaFree(d_list);

create_list(d_list);


print_device_array(d_list,5);


// int *d_list;

// int n = 5;
// cudaMalloc((void**)&d_list,n*sizeof(int));
// int h = 0;
// cudaMemset(d_list,h,n*sizeof(int));
// print_device_array(d_list,n);

cudaDeviceSynchronize();

cudaError_t error = cudaGetLastError();
if(error != cudaSuccess) {
	printf("CUDA error: %s\n",cudaGetErrorString(error));
	exit(-1);
}

}