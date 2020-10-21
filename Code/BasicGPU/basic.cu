typedef struct {
	int x,y;
} INT2;


// __global__ void assign_cell(float2 *dev_cell_list, float2 *dev_unc_pos, int nbead, double lcell, double nbox) {
// 	int global_index = blockIdx.x;

// 	if(global_index < nbead) {
// 		dev_cell_list[global_index].x = floorf(nbox*(dev_unc_pos[global_index].x)/lcell);
// 		dev_cell_list[global_index].y = floorf(nbox*(dev_unc_pos[global_index].y)/lcell);
// 	}
// }
template <class Q>
void arracychck(Q *array, int n) {
Q *temparray = new Q [n];
cudaMemcpy(temparray,array,n*sizeof(Q),cudaMemcpyDeviceToHost);
for(int i = 0 ; i < n ; i++) {
	if(temparray[i]!=temparray[i]) {
		cout << i << endl;
		error("nan array");
	}
}
cout << endl;
delete temparray;

}

template <class Q>
void print_device_array(Q *array, int n) {
Q *temparray = new Q [n];
cudaMemcpy(temparray,array,n*sizeof(Q),cudaMemcpyDeviceToHost);
for(int i = 0 ; i < n ; i++) {
	cout <<setw(4) << temparray[i] <<  ",";
}
cout << endl;
delete temparray;

}


template <class Q>
void print_device_array_indices(Q *array, int n) {
Q *temparray = new Q [n];
cudaMemcpy(temparray,array,n*sizeof(Q),cudaMemcpyDeviceToHost);
for(int i = 0 ; i < n ; i++) {
	if(temparray[i]>0) {
	cout <<setw(4) << i <<  ",";
	}
}
cout << endl;
delete temparray;

}


template <class Q>
void print_device_array(Q *array, int n, int m) {
Q *temparray = new Q [n];
cudaMemcpy(temparray,array,n*sizeof(Q),cudaMemcpyDeviceToHost);
for(int i = 0 ; i < m ; i++) {
	cout << setw(4) << temparray[i] <<  ",";
}
cout << endl;
delete temparray;

}


void print_device_float2(float2 *array, int n) {
float2 *temparray = new float2 [n];
cudaMemcpy(temparray,array,n*sizeof(float2),cudaMemcpyDeviceToHost);
for(int i = 0 ; i < n ; i++) {
	cout << temparray[i].x << " " << temparray[i].y <<  endl;
}
cout << endl;
delete temparray;

}

void file_print_device_float2(float2 *array, int n, ofstream &s) {
float2 *temparray = new float2 [n];
cudaMemcpy(temparray,array,n*sizeof(float2),cudaMemcpyDeviceToHost);
for(int i = 0 ; i < n ; i++) {
	s << temparray[i].x << "," << temparray[i].y <<  endl;
}
delete temparray;

}


void print_device_weave_float2(int *array1, int *array2, double *f1, double *f2, float2 *parts, int n, int N) {
float2 *temparray = new float2 [N];
int *temparray1 = new int [n];
int *temparray2 = new int [n];
double *temparray3 = new double [n];
double *temparray4 = new double [n];

cudaMemcpy(temparray,parts,N*sizeof(float2),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray1,array1,n*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray2,array2,n*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray3,f1,n*sizeof(double),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray4,f2,n*sizeof(double),cudaMemcpyDeviceToHost);

for(int i = 0 ; i < n ; i++) {
	int indx1 = temparray1[i];
	int indx2 = temparray2[i];
	double f1 = temparray3[i];
	double f2 = temparray4[i];

	cout <<"p1: " << indx1 << ": "<<temparray[indx1].x << " " << temparray[indx1].y;
	cout <<" p2: " << indx2 << ": "<<temparray[indx2].x << " " << temparray[indx2].y;
	cout << " :forces: " << f1 << " " << f2 << endl;	
}
cout << endl;
delete temparray;
delete temparray1;
delete temparray2;
delete temparray3;
delete temparray4;
}

void print_device_weave_float2(int *array1, int *array2, double *f1, double *f2, double *f3, double *f4, float2 *parts, int n, int N) {
float2 *temparray = new float2 [N];
int *temparray1 = new int [n];
int *temparray2 = new int [n];
double *temparray3 = new double [n];
double *temparray4 = new double [n];
double *temparray5 = new double [n];
double *temparray6 = new double [n];

cudaMemcpy(temparray,parts,N*sizeof(float2),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray1,array1,n*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray2,array2,n*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray3,f1,n*sizeof(double),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray4,f2,n*sizeof(double),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray5,f3,n*sizeof(double),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray6,f4,n*sizeof(double),cudaMemcpyDeviceToHost);


for(int i = 0 ; i < n ; i++) {
	int indx1 = temparray1[i];
	int indx2 = temparray2[i];
	double f1 = temparray3[i];
	double f2 = temparray4[i];
	double f3 = temparray5[i];
	double f4 = temparray6[i];

	cout <<"p1: " << indx1 << ": "<<temparray[indx1].x << " " << temparray[indx1].y;
	cout <<" p2: " << indx2 << ": "<<temparray[indx2].x << " " << temparray[indx2].y;
	cout << " :forces: " << f1 << " " << f2 << " " << f3 << " " << f4 << endl;	
}
cout << endl;
delete temparray;
delete temparray1;
delete temparray2;
delete temparray3;
delete temparray4;
delete temparray5;
delete temparray6;
}


template <class Q>
void print_device_array_weave(Q *array1, Q *array2, int n) {
Q *temparray1 = new Q [n];
Q *temparray2 = new Q [n];
cudaMemcpy(temparray1,array1,n*sizeof(Q),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray2,array2,n*sizeof(Q),cudaMemcpyDeviceToHost);
for(int i = 0 ; i < n ; i++) {
	cout << temparray1[i] << "," << temparray2[i]<<  "\n";
}
cout << endl;
delete temparray1;
delete temparray2;

}

template <class Q>
void print_device_array_weave(Q *array1, Q *array2, Q *array3, Q * array4, int n) {
Q *temparray1 = new Q [n];
Q *temparray2 = new Q [n];
Q *temparray3 = new Q [n];
Q *temparray4 = new Q [n];

cudaMemcpy(temparray1,array1,n*sizeof(Q),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray2,array2,n*sizeof(Q),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray3,array3,n*sizeof(Q),cudaMemcpyDeviceToHost);
cudaMemcpy(temparray4,array4,n*sizeof(Q),cudaMemcpyDeviceToHost);

for(int i = 0 ; i < n ; i++) {
	cout << temparray1[i] << "," << temparray2[i]<< ","<< temparray3[i]<<"," <<temparray4[i]<<"\n";
}
cout << endl;
delete temparray1;
delete temparray2;
delete temparray3;
delete temparray4;

}

template <class Q>
__global__ void setstate(Q *a, Q i, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		a[global_index]=i;
	}
}

__global__ void setstateincr(int *a, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		int j = global_index;
		//printf("global_index is %d\n", global_index);
		a[global_index]=j;
	}
}
__global__ void setstateincr(int *a, int n, int offset) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		int j = global_index;
		//printf("global_index is %d\n", global_index);
		a[global_index]=j+offset;
	}
}

__global__ void setstate(int *a, int *b, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		a[global_index]=b[global_index];
	}	
}

__global__ void setstaterandomGPU(double *a, double lim, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {

		curandState state;

		curand_init(clock64(),1,0,&state);

		double rr1 = curand_uniform(&state);

		a[global_index] = 2*lim*rr1-lim;		
	}	
}

void setstaterandom(double *a, double lim, int n) {
	setstaterandomGPU<<<n,1>>>(a,lim,n);
}

__global__ void setstatefrom(int *a, int *b,int start, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n&& global_index >= start) {
		a[global_index]=b[global_index-start];
	}
	else if(global_index < start) {
		a[global_index]=0;
	}
}

template <class Q>
__global__ void setstatefunc(int *a, Q func, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		a[global_index]=func(global_index);
	}	
}

template <class Q>
__global__ void setstatefunc(double *a, Q func, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		a[global_index]=func(global_index);
	}	
}

template <class Q>
__global__ void setstatefunc(double *a, int *b,Q func, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		a[global_index]=func(b[global_index]);
	}	
}

__global__ void normalize(double *forcex, double *forcey, double max_s, double v0, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		double mag  = SQR(forcex[global_index])+SQR(forcey[global_index]);
		if(mag > 1E-10) {
		double fac = (mag/(max_s*v0*v0))/tanh(mag/(max_s*v0*v0));
		double rescale = 1./fac;
		forcex[global_index] = rescale * forcex[global_index];
		forcey[global_index] = rescale * forcey[global_index];
		}	
	}
// 	double v0 = (v0_a+v0_b)/2.;
// 	for(int i = 0 ; i < number_of_microtubules ; i++) { //normalize for many motors (not fully collective)
// 		if(motorsbound[i]>1E-10) {
// 			for(int j = 0 ; j < L ; j++) {
// 				double mag =0.0;
// 				for(int k = 0 ; k < dimension ; k++) {
// 					mag+=SQR(forces(na+nb+i*L+j,k));
// 				}
// //				cout << "mag: " << mag << endl;
// 				double fac = (mag/(max_s*v0*v0))/tanh(mag/(max_s*v0*v0));
// //				cout << "fac: " << fac << endl;
				
// 				double rescale = 1./fac;
// 				// cout << "rescaled: " << rescale << endl;
// 				// pausel();
// 				for(int k = 0 ; k < dimension ; k++) {
// 					forces(na+nb+i*L+j,k)=rescale*forces(na+nb+i*L+j,k);
// 				}
// 			}
// 		}
// 	}

}