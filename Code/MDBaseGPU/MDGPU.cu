// typedef struct {
// 	double x,y;
// } float2;
#include "potentialGPU.cu"
#include "potential3GPU.cu"
#include "../BasicGPU/basic.cu"


__device__ void distance_vector2D(float2 h1,float2 h2,float2 &uv, double &d, double l, bool periodic) {
	double dx = h1.x-h2.x;
	double dy = h1.y-h2.y; 
	if(periodic) { 
		double atemp = dx > 0 ? 1. : -1.;
		if(dx*dx > 0.25*l*l) {
			dx = dx - atemp*l;
		}
		double btemp = dy > 0 ? 1. : -1.;
		if(dy*dy > 0.25*l*l) {
			dy = dy - btemp*l;
		}
		// if(dx*dx > l*l ) dx = dx - SIGN(l,dx);
		// if(dy*dy > l*l ) dy = dy - SIGN(l,dy);
	}
	uv.x = dx;
	uv.y = dy;
	d = sqrt(SQR(dx)+SQR(dy));

}

__device__ void correct_position2D(float2 &h1, double l, bool periodic) {
			if(periodic) {
				if(h1.x < 0 ) h1.x = h1.x + l;
				else if(h1.x > l) h1.x = h1.x - l;
				else {
				}
				if(h1.y < 0 ) h1.y = h1.y + l;
				else if(h1.y > l) h1.y = h1.y - l;
				else {
				}				
			}
			else {
				if(h1.x<0) {
					h1.x = -h1.x;
				}
				else if(h1.y>l) {
					h1.y = l-(h1.y-l);
				}
				else{

				}				
			}
}


__global__ void assign_cell(int *dev_cell_list, float2 *dev_unc_pos, int nbead, double lcell, double nbox) {
	//particle positions in dev_unc_pos
	//dev cell list is the list of all particles with their cell
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;

	if(global_index < nbead) {
		dev_cell_list[global_index] = nbox*floorf(nbox*(dev_unc_pos[global_index].x)/lcell) +  floorf(nbox*(dev_unc_pos[global_index].y)/lcell);
		//dev_cell_list[global_index].y =
	}
}

__global__ void cell_counts(int *cell_id, int *cell_counts, int nbead, int cellmax) {
	int global_index = blockIdx.x;
	if(global_index>0 && global_index < nbead ) {
		if(cell_id[global_index]!=cell_id[global_index-1] && cell_id[global_index]==cellmax ) {
			//cell_counts[cell_id[global_index-1]]=global_index;
			//global_index is the point at which the mismatch occurs;
			cell_counts[cell_id[global_index-1]] = global_index;
			cell_counts[cell_id[global_index]] = nbead;
		}
		else if(cell_id[global_index]!=cell_id[global_index-1]){
			cell_counts[cell_id[global_index-1]] = global_index;
		}
		else{

		}
	}
}

//p1 is the list 
__global__ void neighborlist_number(int *p1, int *p2, int *cell_counts, int *c, int nt) {

	int global_index = threadIdx.x + blockIdx.x * blockDim.x;

	if( global_index < nt ) {
		int box1 = p1[global_index];
		int box2 = p2[global_index];



		if(box1 == box2 ) { //if the boxes are the same
			if(box1 == 0 ) { //if it is the first box
				c[global_index] = cell_counts[box1]*(cell_counts[box1]-1)/2;
			}
			else{
				if(cell_counts[box1]==0) { //empty box
					c[global_index]=0;
				}
				else {
					int index1 = cell_counts[box1];
					int box_b = box1 - 1;
					while(box_b >= 0 && cell_counts[box_b] ==0 ) {
						box_b--;
					}
					int pn;
					if(box_b==-1) pn = index1;
					else {

						int index2 = cell_counts[box_b];
						pn = index1 - index2;
					}
					c[global_index] = pn*(pn-1)/2;

				}
			}

		}
		else {
			//box1 is not box2
			if(cell_counts[box1]==0||cell_counts[box2]==0) {
				c[global_index]=0;
			}
			else {
				int index1 = cell_counts[box1];
				int box_b = box1 - 1;
				while(box_b >= 0 && cell_counts[box_b] ==0 ) {
					box_b--;
				}
				int pn1;
				if(box_b==-1) pn1 = index1;
				else {

					int index2 = cell_counts[box_b];
					pn1 = index1 - index2;
				}

				int index3 = cell_counts[box2];
				int box_b2 = box2 - 1;
				while(box_b2 >= 0 && cell_counts[box_b2] ==0 ) {
					box_b2--;
				}
				int pn2;
				if(box_b2==-1) pn2 = index3;
				else {

					int index4 = cell_counts[box_b2];
					pn2 = index3 - index4;
				}
				c[global_index] = pn1*pn2;								
			}

		}
	}

}

__global__ void possible_neighborlist(int *p1, int *p2, int *cell_counts, int *indices_sorted, int *c, int nt, int *i1, int *i2, double *dis2, float2 *positions, double ll, bool periodic) { //store the indices of the possible particles interacting in the list i1,i2
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	//printf("Hello from block %d\n",global_index);
	if( global_index < nt ) {
		// int box1 = p1[global_index];
		// int box2 = p2[global_index];
		// if(global_index==0 && c[global_index] ==0 ) {
		// 	printf("this one %d\n",global_index);
		// }
		if(global_index ==0 && c[0] == 0) {
			//printf("THIS ONE %d\n",global_index);
		}
		else if( c[global_index]-c[global_index-1]==0 ) { //don't do anything if there are no pairs
			//printf("no diff is %d\n", global_index);
			//printf("that one %d\n",global_index);
		}
		else { 
			//printf("other one %d\n",global_index);
			int box1 = p1[global_index];
			int box2 = p2[global_index];
			if(box1 == box2) { //BOXES ARE THE SAME
				int index1 = cell_counts[box1];
				int box_b = box1 - 1;
				while(box_b >= 0 && cell_counts[box_b] ==0 ) {
					box_b--;
				}
				int index2;
				if(box_b==-1) index2 = 0;
				else {
					index2 = cell_counts[box_b];
				}
				int iter = 0;
				int start;
				if(global_index==0) start = 0;
				else start = c[global_index-1];

				// int t_particle1 = indices_sorted[index2];
				// int t_particle2 = indices_sorted[index2+1];
				//printf("box1:%d ,box2: %d, index2: %d, part1: %d, part2: %d, start: %d, nt: %d, p1x: %d, p1y: %d, p2x: %d, p2y: %d, ll: %d, per: %d\n",box1,box2,index2,indices_sorted[index2],indices_sorted[index2+1],start,nt,positions[indices_sorted[index2]].x,positions[indices_sorted[index2]].y,positions[indices_sorted[index2+1]].x,positions[indices_sorted[index2+1]].y,ll,periodic);
				// printf("part1 %d, part2 %d, p1x is %d and p2x is %d, p1y is %d and p2y is %d\n", t_particle1, t_particle2, positions[t_particle1].x , positions[t_particle2].x , positions[t_particle1].y ,positions[t_particle2].y);
				//printf("Hello from global_index %d\n",global_index);
				//printf("Hello from global_index %d, value of start is: %d, iter is %d, index2 is %d and index1 is %d\n",global_index,start,iter,index2,index1);
				//printf("box is %d: ",box1);
				//printf("box: %d box same start is %d, index1 %d, index2 %d\n", box1 , start, index1 ,index2);
				for(int i = index2 ;  i < index1 ; i++) {
					for(int j = i+1 ; j < index1 ; j++) {
							int particle1 =  indices_sorted[i];
							int particle2 =  indices_sorted[j];
							

							//printf("particle1 %d, particle2 %d, start %d and iter %d index1 is %d, j is %d, p1x is %d and p2x is %d, p1y is %d and p2y is %d\n",particle1,particle2,start,iter,index1,j, positions[particle1].x , positions[particle2].x , positions[particle1].y ,positions[particle2].y);
							
							i1[start+iter]=particle1;
							i2[start+iter]=particle2;						


							double lx = positions[particle1].x-positions[particle2].x;

							//double sqrlx =  lx*lx;
							if(periodic) {
								double a = lx > 0. ? 1. : -1.;

								if(lx*lx > 0.25*ll*ll) {
								lx = lx - a*ll;
								}
							}

							double ly = positions[particle1].y-positions[particle2].y;

							//double sqrlx =  lx*lx;
							if(periodic) {
								double b = ly > 0 ? 1. : -1.;

								if(ly*ly > 0.25*ll*ll) {
									ly = ly - b*ll;
								}			
							}				

							dis2[start+iter] = lx*lx+ly*ly;
														


							iter++;
					}
				}



			}
			else{ //BOXES ARE DIFFERENT
				int index1 = cell_counts[box1];
				int box_b = box1 - 1;
				while(box_b >= 0 && cell_counts[box_b] ==0 ) {
					box_b--;
				}
				int index2;
				if(box_b==-1) index2 = 0;
				else {
					index2 = cell_counts[box_b];
				}		

				int index3 = cell_counts[box2];
				int box_b2 = box2 - 1;
				while(box_b2 >= 0 && cell_counts[box_b2] ==0 ) {
					box_b2--;
				}
				int index4;
				if(box_b2==-1) index4 = 0;
				else {
					index4 = cell_counts[box_b2];
				}
				int iter = 0;
				int start;
				if(global_index==0) start = 0;
				else start = c[global_index-1];

				// int t_particle1 = indices_sorted[index2];
				// int t_particle2 = indices_sorted[index4];
				//printf("box1:%d ,box2: %d, index2: %d, index4: %d, part1: %d, part2: %d, start: %d, nt: %d, p1x: %d, p1y: %d, p2x: %d, p2y: %d, ll: %d, per: %d\n",box1,box2,index2,index4,indices_sorted[index2],indices_sorted[index4],start,nt,positions[indices_sorted[index2]].x,positions[indices_sorted[index2]].y,positions[indices_sorted[index4]].x,positions[indices_sorted[index4]].y,ll,periodic);
				// printf("part1 %d, part2 %d, p1x is %d and p2x is %d, p1y is %d and p2y is %d\n", t_particle1, t_particle2, positions[t_particle1].x , positions[t_particle2].x , positions[t_particle1].y ,positions[t_particle2].y);
				//printf("boxes diff %d and %d : ", box1 , box2);
				//printf("box %d, box %d, box different start is %d, index1 %d, index2 %d, index3 %d, index4 %d\n", box1, box2, start, index1, index2, index3 ,index4);
				//printf("Hello from global_index %d, value of start is: %d, iter is %d, index2 is %d and index1 is %d, index 3 is %d and index4 is %d\n",global_index,start,iter,index2,index1,index3,index4);
				for(int i = index2 ;  i < index1 ; i++) {
					for(int j = index4 ; j < index3 ; j++) {	
							int particle1 =  indices_sorted[i];
							int particle2 =  indices_sorted[j];
							//printf("block index %d, thread index %d, particle1 %d, particle2 %d, start %d and iter %d i is %d, j is %d\n",blockIdx.x,threadIdx.x,particle1,particle2,start,iter,i,j);
							//printf("particle1 %d, particle2 %d, start %d and iter %d index1 is %d, index3 is %d, p1x is %d and p2x is %d, p1y is %d and p2y is %d\n",particle1,particle2,start,iter,index1,index3, positions[particle1].x , positions[particle2].x , positions[particle1].y ,positions[particle2].y);
							
							i1[start+iter]=particle1;
							i2[start+iter]=particle2;

							double lx = positions[particle1].x-positions[particle2].x;

							//double sqrlx =  lx*lx;
							if(periodic) {
								double a = lx > 0 ? 1 : -1;

								if(lx*lx > 0.25*ll*ll) {
									lx = lx - a*ll;
								}
							}
							double ly = positions[particle1].y-positions[particle2].y;

							//double sqrlx =  lx*lx;
							if(periodic) {
								double b = ly > 0 ? 1. : -1.;

								if(ly*ly > 0.25*ll*ll) {
									ly = ly - b*ll;
								}							
							}

							dis2[start+iter] = lx*lx+ly*ly;


							iter++;
					}
				}						

			}
			// for(int i = c[global_index-1] ;  i < c[global_index] ; i++) {
			// 	i1[i] = 1;
			// 	i2[i] = 1;

			// }
		}

	}

}

void construct_possible_pair_list(float2 *d_particles, int *d_p_indices, int N, double lcell, int *d_cells1, int *d_cells2, double nbox, bool periodic, int *&d_indices1, int *&d_indices2, double *&d_close, int &tpp2, bool show = false) { //Construct a full pair list with only positions as input
	//assign_cell(int *dev_cell_list, float2 *dev_unc_pos, int nbead, double lcell, double nbox)
	
	//d_particles is the storage of all the particles;
	//d_p_indices is list of all the indices;
	//N is the number of the particles
	//lcell is the length of the 
	//d_cells1,d_cells2 is the list of all the boxes that interact with each other
	//nbox is the number of boxes per 
	//periodic is if there are periodic boundary conditions
	//cut_off is the cut off of the distances

	int st = 5*nbox*nbox;


	int *d_cells;

	int size2 = N*sizeof(int);


	cudaMalloc((void**)&d_cells,size2);

	assign_cell<<<N,1>>>(d_cells,d_particles,N,lcell,nbox); //assign cells to each of the particles in particles, where l is the length of the box and nbox is the number of the boxes per dimension

	//if(show) cout << "got to here 1" << endl;


	thrust::device_ptr<int> t_cells(d_cells);
	thrust::device_ptr<int> t_indices(d_p_indices); //initiliaze device pointers


	// if(show) {
	// cout << "got to here 2" << endl;
	// print_device_array(d_p_indices,N);
	// print_device_array(d_cells,N);
	// cout << N << endl;
	// pausel();
	// }

	// cudaDeviceSynchronize();

	// cudaError_t error = cudaGetLastError();

	// if( error != cudaSuccess) {
	// 	printf("Cuda Error: %s\n",cudaGetErrorString(error));
	// 	exit(-1);
	// }

	thrust::sort_by_key(t_cells,t_cells+N,t_indices); //sort the indices by cell


	//if(show) cout << "got to here 3" << endl;

	thrust::device_vector<int>::iterator iter = thrust::max_element(t_cells,t_cells+N);	

	int largest = *iter;

	int ncells = nbox*nbox;

	int size3 = ncells*sizeof(int);
	int *d_cellsc;
	cudaMalloc((void**)&d_cellsc,size3);
	cudaMemset(d_cellsc,0,size3);


	




	cell_counts<<<N,1>>>(d_cells,d_cellsc,N,largest); //number of particles in each cell

	


	int *d_npb;

	int nbpairs = 5*ncells;
	int size4 = nbpairs*sizeof(int);
	cudaMalloc((void**)&d_npb,size4);
	cudaMemset(d_npb,0,size4);


	neighborlist_number<<<nbpairs,1>>>(d_cells1, d_cells2, d_cellsc,d_npb,nbpairs); //total number of possible pairs




	thrust::device_ptr<int> t_npb(d_npb);
//int tp =  thrust::reduce(t_npb,t_npb+nbpairs);

	thrust::inclusive_scan(t_npb,t_npb+nbpairs,t_npb); //cumulative binnings




	int tpp;
	cudaMemcpy(&tpp,d_npb+nbpairs-1,sizeof(int),cudaMemcpyDeviceToHost);
	//cout << tpp << endl;
	tpp2 = tpp;
	// int tpp = *d_tpp;
	int size5 = tpp*sizeof(int);





	cudaMalloc((void**)&d_indices1,tpp*sizeof(int));

	cudaMalloc((void**)&d_indices2,tpp*sizeof(int));

	cudaMalloc((void**)&d_close,tpp*sizeof(double));

	setstate<<<tpp,1>>>(d_indices1,0,tpp);

	setstate<<<tpp,1>>>(d_indices2,0,tpp);

	setstate<<<tpp,1>>>(d_close,1.5,tpp);

	possible_neighborlist<<<nbpairs,1>>>(d_cells1, d_cells2, d_cellsc, d_p_indices, d_npb, nbpairs, d_indices1, d_indices2, d_close, d_particles,lcell,periodic);


	if(show) cout << "device delete" << endl;
	cudaFree(d_cellsc); 
	if(show) cout << "cuda free 1" << endl;
	cudaFree(d_npb); 
	if(show) cout << "cuda free 2" << endl;
	cudaFree(d_cells);
	if(show) cout << "cuda free 3" << endl;



}

struct less_than_condition {
double dis_less;
__host__ __device__ less_than_condition(double dis) : dis_less(dis) {}
__host__ __device__ int operator()(int i1, int i2, double a) { 
	if(a < dis_less) return 1;
	else return 0;
}
};



struct less_than_condition_NAND {
	double dis_less;
	int p1,p2;
	int p3,p4;
	__host__ __device__ less_than_condition_NAND(double dis, int q1, int q2, int q3, int q4) : dis_less(dis),p1(q1),p2(q2),p3(q3),p4(q4) {}
	__host__ __device__ int operator()(int i1, int i2, double a) { 
		int j1,j2;
		if(i1 < i2) {j1 = i1; j2 = i2;}
		else{j1 = i2 ; j2 = i1; }

		if(j1<p2&&j1>=p1) {
			if(j2<p4&&j2>=p3){
				if(a < dis_less) {
					return 1;
				}
				else{
					return 0;
				}
			}
			else{
				return 0; 
			}
		}
		else{
			return 0;
		}

	}
};

struct less_than_condition_AND {
	double dis_less;
	int p1,p2;
	__host__ __device__ less_than_condition_AND(double dis, int q1, int q2) : dis_less(dis),p1(q1),p2(q2) {}
	__host__ __device__ int operator()(int i1, int i2, double a) {
		if(i1<p2&&i1>=p1) {
			if(i2<p2&&i2>=p1){
				if(a < dis_less) {
					return 1;
				}
				else{
					return 0;
				}
			}
			else{
				return 0; 
			}
		}
		else{
			return 0;
		}
	}
};

template <typename cond>
__global__ void applycondition(int *d_indices1, int *d_indices2,double *d_close,cond F,int N, int *res) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < N ) {
		int i1 = d_indices1[global_index];
		int i2 = d_indices2[global_index];
		double dis =  d_close[global_index];

		res[global_index]=F(i1,i2,dis);
	}
}

__global__ void neighborlist(int *p1, int *p2, int *dis, int nt, int *i1, int *i2, int *i3, int *i4) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;

	if(global_index ==0 ) {
		if(dis[global_index]==0) {

		}
		else{
			i1[0] = p1[global_index];
			i2[0] = p2[global_index];
			i3[0] = p1[global_index];
			i4[0] = p2[global_index];
		}

	}
	else if(global_index < nt ) {
		if(dis[global_index]!=dis[global_index-1]) {
			int fg = dis[global_index-1];
			i1[fg] = p1[global_index];
			i2[fg] = p2[global_index];
			i3[fg] = p1[global_index];
			i4[fg] = p2[global_index];			

		}

	}
	else {

	}
}


template <typename Func>
void pairlist(int *d_indices1,int *d_indices2, double *d_close,Func f,int *&d_list1,int *&d_list2,int *&d_list3,int *&d_list4, int tpp, int &th2, int show = true) {

int *d_close2;
cudaMalloc((void**)&d_close2,tpp*sizeof(int));
//cudaMemset(d_close2,0,tpp*sizeof(int));


setstate<<<tpp,1>>>(d_close2,0,tpp);

//resetstate<<<tpp,1>>>(d_close2,0,tpp);

//if(show) cout << "mem allocated" << endl;

applycondition<<<tpp,1>>>(d_indices1,d_indices2,d_close,f,tpp,d_close2);

//if(show) cout << "cond applied" << endl;

//if(show) print_device_array(d_close,tpp);


thrust::device_ptr<int> t_close(d_close2);

thrust::inclusive_scan(t_close,t_close+tpp,t_close);

//if(show) cout << "condscan" << endl;

int th;

cudaMemcpy(&th,d_close2+tpp-1,sizeof(int),cudaMemcpyDeviceToHost);

//if(show) cout << "memcpy" << endl;

th2 = th;

cudaMalloc((void**)&d_list1,th*sizeof(int));

cudaMalloc((void**)&d_list2,th*sizeof(int));

cudaMalloc((void**)&d_list3,th*sizeof(int));

cudaMalloc((void**)&d_list4,th*sizeof(int));

neighborlist<<<tpp,1>>>(d_indices1, d_indices2, d_close2, tpp, d_list1, d_list2,d_list3,d_list4);

cudaFree(d_close2);


}


// __global__ void neighborlist(bool *list, int *p1, int *p2) {
// 	//p1 and p2 are a list of indices of particles which are in the correct part
// }
//given a list of 


template <typename Q>
__global__ void calculateforces2DGPU(int *i1, int *i2, float2 *positions, double *forces1x, double *forces1y, double *forces2x, double *forces2y, Q iny ,int nt, double ll, bool periodic) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;

	if(global_index < nt) {
		int particle1 = i1[global_index];
		int particle2 = i2[global_index];

		double lx = positions[particle1].x-positions[particle2].x;

		//double sqrlx =  lx*lx;
		if(periodic) {
			double a = lx > 0 ? 1 : -1;
			if(lx*lx > 0.25*ll*ll) {
				lx = lx - a*ll;
			}
		}	

		double ly = positions[particle1].y-positions[particle2].y;

		//double sqrlx =  lx*lx;
		if(periodic) {
			double b = ly > 0 ? 1 : -1;
			if(ly*ly > 0.25*ll*ll) {
				ly = ly - b*ll;
			}
		}


		double dis = sqrt(lx*lx+ly*ly);

		double f = iny(dis);

		forces1x[global_index] = lx*f/dis;
		forces1y[global_index] = ly*f/dis;

		forces2x[global_index] = -lx*f/dis;
		forces2y[global_index] = -ly*f/dis;

	}


}


template <typename Q>
__global__ void calculateforces_threebodyGPU(float2 *particles, int *list1, int *list2, int *list3,double *forcep1x, double *forcep1y, double *forcep2x, double *forcep2y,double *forcep3x, double *forcep3y,Q iny,double ll, bool periodic, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < n) {
		int p1 = list1[global_index];
		int p2 = list2[global_index];
		int p3 = list3[global_index];


		float2 a = particles[p1];

		float2 b = particles[p2];

		float2 c = particles[p3];

		float2 ab;

		float2 bc;

		double lx = b.x-a.x;

		//double sqrlx =  lx*lx;
		if(periodic) {
			double atemp = lx > 0. ? 1. : -1.;

			if(lx*lx > 0.25*ll*ll) {
				lx = lx - atemp*ll;
			}
		}
		double ly = b.y-a.y;

		//double sqrlx =  lx*lx;
		if(periodic) {
			double btemp = ly > 0. ? 1. : -1.;

			if(ly*ly > 0.25*ll*ll) {
				ly = ly - btemp*ll;
			}							
		}		
		ab.x = lx;
		ab.y = ly;

		double lx2 = c.x-b.x;

		//double sqrlx =  lx*lx;
		if(periodic) {
			double atemp = lx2 > 0. ? 1. : -1.;

			if(lx2*lx2 > 0.25*ll*ll) {
				lx2 = lx2 - atemp*ll;
			}
		}
		double ly2 = c.y-b.y;

		//double sqrlx =  lx*lx;
		if(periodic) {
			double btemp = ly2 > 0 ? 1 : -1;

			if(ly2*ly2 > 0.25*ll*ll) {
				ly2 = ly2 - btemp*ll;
			}							
		}		
		bc.x = lx2;
		bc.y = ly2;		


		float2_3 forc = iny(ab,bc);

		//printf("gin: %d, 1: %f 2: %f 3: %f 4: %f 5: %f 6: %f, 7: %f, 8: %f, 9: %f, 10: %f\n",global_index,forc.f1.x,forc.f1.y,forc.f2.x,forc.f2.y,forc.f3.x,forc.f3.y,ab.x,ab.y,bc.x,bc.y);
		//float2_3 force(float2 ab,float2 bc)

		forcep1x[global_index] = forc.f1.x;

		forcep1y[global_index] = forc.f1.y;

		forcep2x[global_index] = forc.f2.x;

		forcep2y[global_index] = forc.f2.y;

		forcep3x[global_index] = forc.f3.x;

		forcep3y[global_index] = forc.f3.y;




	}

}


template <typename Q>
void calculateforces2D(int *i1, int *i2, float2 *positions, double *&forces1x, double *&forces1y, double *&forces2x, double *&forces2y, Q iny ,int nt, double ll, bool periodic) {

calculateforces2DGPU<<<nt,1>>>(i1,i2,positions,forces1x,forces1y,forces2x,forces2y,iny ,nt,ll,periodic);

}

__global__ void addforce(double *forcex, double *forcey, int *i1, int *i2, double *forces1x, double *forces1y, int nt) {
	 int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	 if(global_index < nt) {
	 	int particle1 = i1[global_index];
		int particle2 = i2[global_index];

		forcex[particle1] += forces1x[global_index];
		forcey[particle2] += forces1y[global_index];
	 }
}


//take the forces and sum them

void ReduceForces(int *dd_list1,int *dd_list2,int *dd_list3,int *dd_list4,double *d_forces1x,double *d_forces2x,double *d_forces1y,double *d_forces2y,double *d_totalforcex,double *d_totalforcey, int th) {

int *d_list1;
int *d_list2;
int *d_list3;
int *d_list4;
cudaMalloc((void**)&d_list1,th*sizeof(double));
cudaMalloc((void**)&d_list2,th*sizeof(double));
cudaMalloc((void**)&d_list3,th*sizeof(double));
cudaMalloc((void**)&d_list4,th*sizeof(double));

setstate<<<th,1>>>(d_list1, dd_list1, th);
setstate<<<th,1>>>(d_list2, dd_list2, th);
setstate<<<th,1>>>(d_list3, dd_list3, th);
setstate<<<th,1>>>(d_list4, dd_list4, th);



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

new_end1 = thrust::reduce_by_key(t_list1,t_list1+th,t_forces1x,t_key_reduce1x,t_sumforces1x);


new_end2 = thrust::reduce_by_key(t_list3,t_list3+th,t_forces1y,t_key_reduce1y,t_sumforces1y);


new_end3 = thrust::reduce_by_key(t_list2,t_list2+th,t_forces2x,t_key_reduce2x,t_sumforces2x);


new_end4 = thrust::reduce_by_key(t_list4,t_list4+th,t_forces2y,t_key_reduce2y,t_sumforces2y);




int ih1 = thrust::raw_pointer_cast(&new_end1.first[0])-thrust::raw_pointer_cast(&t_key_reduce1x[0]);

int ih3 = thrust::raw_pointer_cast(&new_end3.first[0])-thrust::raw_pointer_cast(&t_key_reduce2x[0]);

addforce<<<ih1,1>>>(d_totalforcex,d_totalforcey, d_key_reduce1x,d_key_reduce1y, d_sumforces1x, d_sumforces1y, ih1);

addforce<<<ih3,1>>>(d_totalforcex,d_totalforcey, d_key_reduce2x,d_key_reduce2y, d_sumforces2x, d_sumforces2y, ih3);

cudaFree(d_sumforces1x);
cudaFree(d_sumforces2x);
cudaFree(d_sumforces1y);
cudaFree(d_sumforces2y);


cudaFree(d_list1);
cudaFree(d_list2);
cudaFree(d_list3);
cudaFree(d_list4);

cudaFree(d_key_reduce1x);
cudaFree(d_key_reduce2x);
cudaFree(d_key_reduce1y);
cudaFree(d_key_reduce2y);

}

void ReduceForces(int *dd_list1,int *dd_list2,double *d_forces1x,double *d_forces2x,double *d_forces1y,double *d_forces2y,double *d_totalforcex,double *d_totalforcey, int th) {

int *d_list1;
int *d_list2;
int *d_list3;
int *d_list4;
cudaMalloc((void**)&d_list1,th*sizeof(double));
cudaMalloc((void**)&d_list2,th*sizeof(double));
cudaMalloc((void**)&d_list3,th*sizeof(double));
cudaMalloc((void**)&d_list4,th*sizeof(double));

setstate<<<th,1>>>(d_list1, dd_list1, th);
setstate<<<th,1>>>(d_list2, dd_list2, th);
setstate<<<th,1>>>(d_list3, dd_list1, th);
setstate<<<th,1>>>(d_list4, dd_list2, th);



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

new_end1 = thrust::reduce_by_key(t_list1,t_list1+th,t_forces1x,t_key_reduce1x,t_sumforces1x);


new_end2 = thrust::reduce_by_key(t_list3,t_list3+th,t_forces1y,t_key_reduce1y,t_sumforces1y);


new_end3 = thrust::reduce_by_key(t_list2,t_list2+th,t_forces2x,t_key_reduce2x,t_sumforces2x);


new_end4 = thrust::reduce_by_key(t_list4,t_list4+th,t_forces2y,t_key_reduce2y,t_sumforces2y);




int ih1 = thrust::raw_pointer_cast(&new_end1.first[0])-thrust::raw_pointer_cast(&t_key_reduce1x[0]);

int ih3 = thrust::raw_pointer_cast(&new_end3.first[0])-thrust::raw_pointer_cast(&t_key_reduce2x[0]);

addforce<<<ih1,1>>>(d_totalforcex,d_totalforcey, d_key_reduce1x,d_key_reduce1y, d_sumforces1x, d_sumforces1y, ih1);

addforce<<<ih3,1>>>(d_totalforcex,d_totalforcey, d_key_reduce2x,d_key_reduce2y, d_sumforces2x, d_sumforces2y, ih3);

cudaFree(d_sumforces1x);
cudaFree(d_sumforces2x);
cudaFree(d_sumforces1y);
cudaFree(d_sumforces2y);

cudaFree(d_list1);
cudaFree(d_list2);
cudaFree(d_list3);
cudaFree(d_list4);

cudaFree(d_key_reduce1x);
cudaFree(d_key_reduce2x);
cudaFree(d_key_reduce1y);
cudaFree(d_key_reduce2y);

}


void ReduceForcesAndNormalize(int *dd_list1,double *d_forces1x,double *d_forces1y,double *d_totalforcex,double *d_totalforcey, double max_s, double v0, int th) {

int *d_list1;
int *d_list2;
cudaMalloc((void**)&d_list1,th*sizeof(double));
cudaMalloc((void**)&d_list2,th*sizeof(double));

setstate<<<th,1>>>(d_list1, dd_list1, th);
setstate<<<th,1>>>(d_list2, dd_list1, th);



thrust::device_ptr<double> t_forces1x(d_forces1x);
thrust::device_ptr<double> t_forces1y(d_forces1y);
thrust::device_ptr<int> t_list1(d_list1);
thrust::device_ptr<int> t_list2(d_list2);

// thrust::device_ptr<double> t_forces2x(d_forces2x);
// thrust::device_ptr<double> t_forces2y(d_forces2y);
// thrust::device_ptr<int> t_list2(d_list2);
// thrust::device_ptr<int> t_list4(d_list4);

thrust::sort_by_key(t_list1,t_list1+th,t_forces1x);


thrust::sort_by_key(t_list2,t_list2+th,t_forces1y);


double *d_sumforces1x;
double *d_sumforces1y;

int *d_key_reduce1x;
int *d_key_reduce1y;

cudaMalloc((void**)&d_sumforces1x,th*sizeof(double));
cudaMalloc((void**)&d_sumforces1y,th*sizeof(double));


cudaMalloc((void**)&d_key_reduce1x,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce1y,th*sizeof(int));


thrust::device_ptr<double> t_sumforces1x(d_sumforces1x);
thrust::device_ptr<double> t_sumforces1y(d_sumforces1y);

thrust::device_ptr<int> t_key_reduce1x(d_key_reduce1x);
thrust::device_ptr<int> t_key_reduce1y(d_key_reduce1y);



thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end1;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end2;


new_end1 = thrust::reduce_by_key(t_list1,t_list1+th,t_forces1x,t_key_reduce1x,t_sumforces1x);


new_end2 = thrust::reduce_by_key(t_list2,t_list2+th,t_forces1y,t_key_reduce1y,t_sumforces1y);





int ih1 = thrust::raw_pointer_cast(&new_end1.first[0])-thrust::raw_pointer_cast(&t_key_reduce1x[0]);



normalize<<<ih1,1>>>(d_sumforces1x, d_sumforces1y,max_s, v0, ih1);


addforce<<<ih1,1>>>(d_totalforcex,d_totalforcey, d_key_reduce1x,d_key_reduce1y, d_sumforces1x, d_sumforces1y, ih1);


cudaFree(d_sumforces1x);
cudaFree(d_sumforces1y);

cudaFree(d_list1);
cudaFree(d_list2);

cudaFree(d_key_reduce1x);
cudaFree(d_key_reduce1y);


}



void ReduceForces3(int *dd_list1,int *dd_list2,int *dd_list3,double *d_forces1x,double *d_forces2x, double *d_forces3x,double *d_forces1y,double *d_forces2y, double *d_forces3y,double *d_totalforcex,double *d_totalforcey, int th) {

int *d_list1;
int *d_list2;
int *d_list3;
int *d_list4;
int *d_list5;
int *d_list6;
cudaMalloc((void**)&d_list1,th*sizeof(double));
cudaMalloc((void**)&d_list2,th*sizeof(double));
cudaMalloc((void**)&d_list3,th*sizeof(double));
cudaMalloc((void**)&d_list4,th*sizeof(double));
cudaMalloc((void**)&d_list5,th*sizeof(double));
cudaMalloc((void**)&d_list6,th*sizeof(double));



setstate<<<th,1>>>(d_list1, dd_list1, th);
setstate<<<th,1>>>(d_list2, dd_list2, th);
setstate<<<th,1>>>(d_list3, dd_list3, th);
setstate<<<th,1>>>(d_list4, dd_list1, th);
setstate<<<th,1>>>(d_list5, dd_list2, th);
setstate<<<th,1>>>(d_list6, dd_list3, th);




thrust::device_ptr<double> t_forces1x(d_forces1x);
thrust::device_ptr<double> t_forces1y(d_forces1y);
thrust::device_ptr<int> t_list1(d_list1);
thrust::device_ptr<int> t_list4(d_list4);

thrust::device_ptr<double> t_forces2x(d_forces2x);
thrust::device_ptr<double> t_forces2y(d_forces2y);
thrust::device_ptr<int> t_list2(d_list2);
thrust::device_ptr<int> t_list5(d_list5);

thrust::device_ptr<double> t_forces3x(d_forces3x);
thrust::device_ptr<double> t_forces3y(d_forces3y);
thrust::device_ptr<int> t_list3(d_list3);
thrust::device_ptr<int> t_list6(d_list6);

thrust::sort_by_key(t_list1,t_list1+th,t_forces1x);


thrust::sort_by_key(t_list2,t_list2+th,t_forces2x);

thrust::sort_by_key(t_list3,t_list3+th,t_forces3x);

thrust::sort_by_key(t_list4,t_list4+th,t_forces1y);

thrust::sort_by_key(t_list5,t_list5+th,t_forces2y);

thrust::sort_by_key(t_list6,t_list6+th,t_forces3y);

double *d_sumforces1x;
double *d_sumforces1y;
double *d_sumforces2x;
double *d_sumforces2y;
double *d_sumforces3x;
double *d_sumforces3y;


int *d_key_reduce1x;
int *d_key_reduce1y;
int *d_key_reduce2x;
int *d_key_reduce2y;
int *d_key_reduce3x;
int *d_key_reduce3y;

cudaMalloc((void**)&d_sumforces1x,th*sizeof(double));
cudaMalloc((void**)&d_sumforces2x,th*sizeof(double));
cudaMalloc((void**)&d_sumforces1y,th*sizeof(double));
cudaMalloc((void**)&d_sumforces2y,th*sizeof(double));
cudaMalloc((void**)&d_sumforces3x,th*sizeof(double));
cudaMalloc((void**)&d_sumforces3y,th*sizeof(double));



cudaMalloc((void**)&d_key_reduce1x,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce2x,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce1y,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce2y,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce3x,th*sizeof(int));
cudaMalloc((void**)&d_key_reduce3y,th*sizeof(int));

thrust::device_ptr<double> t_sumforces1x(d_sumforces1x);
thrust::device_ptr<double> t_sumforces2x(d_sumforces2x);
thrust::device_ptr<double> t_sumforces1y(d_sumforces1y);
thrust::device_ptr<double> t_sumforces2y(d_sumforces2y);
thrust::device_ptr<double> t_sumforces3x(d_sumforces3x);
thrust::device_ptr<double> t_sumforces3y(d_sumforces3y);

thrust::device_ptr<int> t_key_reduce1x(d_key_reduce1x);
thrust::device_ptr<int> t_key_reduce1y(d_key_reduce1y);
thrust::device_ptr<int> t_key_reduce2x(d_key_reduce2x);
thrust::device_ptr<int> t_key_reduce2y(d_key_reduce2y);
thrust::device_ptr<int> t_key_reduce3x(d_key_reduce3x);
thrust::device_ptr<int> t_key_reduce3y(d_key_reduce3y);


thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end1;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end2;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end3;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end4;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end5;

thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<double> > new_end6;


new_end1 = thrust::reduce_by_key(t_list1,t_list1+th,t_forces1x,t_key_reduce1x,t_sumforces1x);


new_end2 = thrust::reduce_by_key(t_list4,t_list4+th,t_forces1y,t_key_reduce1y,t_sumforces1y);


new_end3 = thrust::reduce_by_key(t_list2,t_list2+th,t_forces2x,t_key_reduce2x,t_sumforces2x);


new_end4 = thrust::reduce_by_key(t_list5,t_list5+th,t_forces2y,t_key_reduce2y,t_sumforces2y);


new_end5 = thrust::reduce_by_key(t_list3,t_list3+th,t_forces3x,t_key_reduce3x,t_sumforces3x);


new_end6 = thrust::reduce_by_key(t_list6,t_list6+th,t_forces3y,t_key_reduce3y,t_sumforces3y);




int ih1 = thrust::raw_pointer_cast(&new_end1.first[0])-thrust::raw_pointer_cast(&t_key_reduce1x[0]);

int ih3 = thrust::raw_pointer_cast(&new_end3.first[0])-thrust::raw_pointer_cast(&t_key_reduce2x[0]);

int ih5 = thrust::raw_pointer_cast(&new_end5.first[0])-thrust::raw_pointer_cast(&t_key_reduce3x[0]);

addforce<<<ih1,1>>>(d_totalforcex,d_totalforcey, d_key_reduce1x,d_key_reduce1y, d_sumforces1x, d_sumforces1y, ih1);

addforce<<<ih3,1>>>(d_totalforcex,d_totalforcey, d_key_reduce2x,d_key_reduce2y, d_sumforces2x, d_sumforces2y, ih3);

addforce<<<ih3,1>>>(d_totalforcex,d_totalforcey, d_key_reduce3x,d_key_reduce3y, d_sumforces3x, d_sumforces3y, ih5);

cudaFree(d_sumforces1x);

cudaFree(d_sumforces2x);

cudaFree(d_sumforces1y);

cudaFree(d_sumforces2y);

cudaFree(d_sumforces3x);

cudaFree(d_sumforces3y);

cudaFree(d_list1);
cudaFree(d_list2);
cudaFree(d_list3);
cudaFree(d_list4);
cudaFree(d_list5);
cudaFree(d_list6);

cudaFree(d_key_reduce1x);
cudaFree(d_key_reduce2x);
cudaFree(d_key_reduce1y);
cudaFree(d_key_reduce2y);
cudaFree(d_key_reduce3x);
cudaFree(d_key_reduce3y);


}
  
__global__ void advmom(double *p, double *F, double *R, int nt, double cons1, double cons2, double cons3) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x; 
	if(global_index < nt ) {
		p[global_index] =  cons1*p[global_index] + cons2*F[global_index]+cons3*R[global_index];
	} 
}

__global__ void advpos(double *x, double *p, int nt, double cons1) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x; 
	if(global_index < nt ) {
		//(*dat)(i,i1) = (*dat)(i,i1)+ c1*(*mom)(i,i1);
		x[global_index] =  x[global_index] + cons1*p[global_index];
	} 
}


__global__ void applypbc(double *x, double *p, double l, bool periodic, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x; 
	if(global_index < n) {
		if(periodic) {
			if(x[global_index] < 0) {
				x[global_index] = x[global_index]+l;
			}
			else if(x[global_index]>l) {
				x[global_index] = x[global_index]-l;
			}
			else{

			}
		}
		else{
			if(x[global_index] < 0) {
				x[global_index] = -x[global_index];
				p[global_index] = -p[global_index];
			}
			else if(x[global_index]>l) {
				x[global_index] = l-(x[global_index]-l);
				p[global_index] = -p[global_index];
			}	
			else{

			}	
		}
	}	
}

__global__ void applypbc2DGPU(float2 *x, float2 *p, double l, bool periodic, int n) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x; 
	if(global_index < n) {
		if(periodic) {
			if(x[global_index].x < 0) {
				x[global_index].x = x[global_index].x+l;
			}
			else if(x[global_index].x>l) {
				x[global_index].x = x[global_index].x-l;
			}
			else{

			}
			if(x[global_index].y < 0) {
				x[global_index].y = x[global_index].y+l;
			}
			else if(x[global_index].y>l) {
				x[global_index].y = x[global_index].y-l;
			}
			else{

			}			
		}
		else{
			if(x[global_index].x < 0) {
				x[global_index].x = -x[global_index].x;
				p[global_index].x = -p[global_index].x;
			}
			else if(x[global_index].x>l) {
				x[global_index].x = l-(x[global_index].x-l);
				p[global_index].x = -p[global_index].x;
			}	
			else{

			}	
			if(x[global_index].y < 0) {
				x[global_index].y = -x[global_index].y;
				p[global_index].y = -p[global_index].y;
			}
			else if(x[global_index].y>l) {
				x[global_index].y = l-(x[global_index].y-l);
				p[global_index].y = -p[global_index].y;
			}	
			else{

			}
		}
	}	
}

void applypbc2D(float2 *x, float2 *p, double l, bool periodic, int n) {
applypbc2DGPU<<<n,1>>>(x,p,l,periodic,n);
}


__global__ void advmom2DGPU(float2 *p, double *Fx, double *Fy, double *Rx, double *Ry, double cons1, double cons2, double cons3, int nt) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x; 
	if(global_index < nt ) {
		p[global_index].x =  cons1*p[global_index].x + cons2*Fx[global_index]+cons3*Rx[global_index];
		p[global_index].y =  cons1*p[global_index].y + cons2*Fy[global_index]+cons3*Ry[global_index];		
	} 
}

template <typename Func>
__global__ void advmom2DGPU_spatialdependence(float2 *p, float2 *x, double *Fx, double *Fy, double *Rx, double *Ry, Func fun, double dt, double kT, double m, int nt){ 
	int global_index = threadIdx.x + blockIdx.x * blockDim.x; 
	if(global_index < nt ) {

		double p_x = x[global_index].x;
		double p_y = x[global_index].y;
		double t_gamma = fun(p_x,p_y);
		double t_d = (t_gamma*dt/2.);
		double t_q = (dt)/2.;
		double t_r = sqrt(0.5*kT*(t_gamma)*(m)*(dt));

		double t_c2 = (1.0/(1.0+(t_d)));
		double t_c3 = (1.0/(1.0+(t_d)))*t_q;
		double t_c4 = (1.0/(1.0+(t_d)))*t_r;
		double t_c5 = (1-(t_d));
		p[global_index].x =  (t_c5*t_c2)*p[global_index].x + (t_c5*(t_c3)+t_q)*Fx[global_index]+(t_c5*(t_c4)+t_r)*Rx[global_index];
		p[global_index].y =  (t_c5*t_c2)*p[global_index].y + (t_c5*(t_c3)+t_q)*Fy[global_index]+(t_c5*(t_c4)+t_r)*Ry[global_index];		
	} 	
}

template <typename Func>
__global__ void advmom2DGPU_particledependence(float2 *p, double *Fx, double *Fy, double *Rx, double *Ry, Func fun, double dt, double kT, double m, int nt){ 
	int global_index = threadIdx.x + blockIdx.x * blockDim.x; 
	if(global_index < nt ) {


		double t_gamma = fun(global_index);//fun(p_x,p_y);
		double t_d = (t_gamma*dt/2.);
		double t_q = (dt)/2.;
		double t_r = sqrt(0.5*kT*(t_gamma)*(m)*(dt));

		double t_c2 = (1.0/(1.0+(t_d)));
		double t_c3 = (1.0/(1.0+(t_d)))*t_q;
		double t_c4 = (1.0/(1.0+(t_d)))*t_r;
		double t_c5 = (1-(t_d));
		p[global_index].x =  (t_c5*t_c2)*p[global_index].x + (t_c5*(t_c3)+t_q)*Fx[global_index]+(t_c5*(t_c4)+t_r)*Rx[global_index];
		p[global_index].y =  (t_c5*t_c2)*p[global_index].y + (t_c5*(t_c3)+t_q)*Fy[global_index]+(t_c5*(t_c4)+t_r)*Ry[global_index];		
	} 	
}


__global__ void advpos2DGPU(float2 *x, float2 *p, double cons1, int nt) {
	int global_index = threadIdx.x + blockIdx.x * blockDim.x; 
	if(global_index < nt ) {
		//(*dat)(i,i1) = (*dat)(i,i1)+ c1*(*mom)(i,i1);
		x[global_index].x =  x[global_index].x + cons1*p[global_index].x;
		x[global_index].y =  x[global_index].y + cons1*p[global_index].y;
	} 
}

template <typename Func>
void advmom2D_spatialdependence(float2 *p, float2 *x, double *Fx, double *Fy, double *Rx, double *Ry, Func fun, double dt, double kT, double m, int nt) { 
advmom2DGPU_spatialdependence<<<nt,1>>>(p,x,Fx,Fy,Rx,Ry,fun,dt,kT,m, nt);
}

template <typename Func>
void advmom2D_particledependence(float2 *p, double *Fx, double *Fy, double *Rx, double *Ry, Func fun, double dt, double kT, double m, int nt) { 
advmom2DGPU_particledependence<<<nt,1>>>(p,Fx,Fy,Rx,Ry,fun,dt,kT,m, nt);
}

void advmom2D(float2 *p, double *Fx, double *Fy, double *Rx, double *Ry, double cons1, double cons2, double cons3, int nt) {
advmom2DGPU<<<nt,1>>>(p,Fx,Fy,Rx,Ry,cons1,cons2,cons3, nt);
}

void advpos2D(float2 *x, float2 *p, double cons1, int nt) {
advpos2DGPU<<<nt,1>>>(x,p,cons1,nt);
}