__global__ void CalculateUnbindingsGPU(float2 *particles, int *bound, double *bound_along, double probunbind, int na, int nb, double L, int N, int *changestate, double l, bool periodic, double disless) { //bindings has binding info, 
	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < N) {
		if(bound[global_index]>0) {
			double tt  =  bound_along[global_index];
			int mt =  bound[global_index];
			double tt2 = tt*(double)(L-1);


			int p2 = global_index;
			int k2 = ceil(tt2);
			int k1 = k2 - 1;
			double t = tt2-k1;
			int p1a = na+nb+(mt-1)*L+k1;
			int p1b = na+nb+(mt-1)*L+k2;
			//vector1<double> a0(dimension),a1(dimension),a2(dimension);

			float2 a0,a1,a2;

			a0.x=particles[p2].x;
			a0.y=particles[p2].y;
			a1.x=particles[p1a].x;
			a1.y=particles[p1a].y;
			a2.x=particles[p1b].x;
			a2.y=particles[p1b].y;


			//vector1<double> uv3(dimension),uv4(dimension);
			double dis3,dis4;
			float2 uv3,uv4;

			distance_vector2D(a2,a1,uv3,dis3,l,periodic);

			float2 tempvector;
			tempvector.x = a1.x + t*uv3.x;
			tempvector.y = a1.y + t*uv3.y;

			correct_position2D(tempvector,l,periodic);

			curandState state;

			curand_init(clock64(),1,0,&state);

			double rr1 = curand_uniform(&state);

			distance_vector2D(a0,tempvector,uv4,dis4,l,periodic);

				if(abs(bound_along[global_index]-0.5)>0.4) { //drop off
					changestate[global_index] = 1;
					bound[global_index]=0;
					bound_along[global_index]=0.;
				}
				else if(rr1<probunbind) { //unbind with rate
					changestate[global_index] = 1;
					bound[global_index]=0;
					bound_along[global_index]=0.;
				}
				else if(dis4>1.2*disless) { // if distance is larger
					changestate[global_index] = 1;
					bound[global_index]=0;
					bound_along[global_index]=0.;				
				}
				else {

				} 
			// vector1<double> tempvector = a1+(t*(uv3));
			// (obj->getgeo()).correct_position(tempvector);
			// (obj->getgeo()).distance_vector(a0,tempvector,uv4,dis4);
		}
	}
}

void Microtubule::callCalculateUnbindingsGPU(float2 *particles, int *bound, double *bound_along, int *changestate) {
	int NaNb = na+nb;
	CalculateUnbindingsGPU<<<NaNb,1>>>(particles,bound,bound_along,probunbind,na,nb,L,totalN,changestate,l,is_periodic,excess_force_distance);
}


__device__ void BindGPU(int p1, int p2, float2 h1, float2 h2, int na, int nb, double L, double l, bool periodic, double disless, int initialbind, int &finalbind, double &fpos, double probbind) { //p1 is the microtubule
		int mt = 1+(int)((double)(p1-na-nb)/L); //which microtubule, starting from 1
	//	cout << 1 << endl;
		
		int LL = (int)L;
		double tt = (double)((p1-na-nb)%LL);
	//	cout << 2 << endl;
		tt = tt/(double)(L-1);
	//	cout << 3 << endl;
		// double dx = h1.x-h2.x;
		// double dy = h1.y-h2.y; 
		// if(periodic) { 
		// 	if(SQR(dx) > l*l ) dx = dx - SIGN(l,dx);
		// 	if(SQR(dy) > l*l ) dy = dy - SIGN(l,dy);
		// }
		// double dis4 = sqrt(SQR(dx)+SQR(dy));		
		float2 uv;
		double dis4;
		distance_vector2D(h1,h2,uv, dis4,l,periodic);
	//	cout << 4 << endl;
		//printf("p1: %d, p2: %d, initialbind: %d, (x1,y1) %d,%d, (x2,y2) %d,%d, tt %d\n",p1,p2,initialbind,h1.x,h1.y,h2.x,h2.y,tt);

		if(initialbind==0 && abs(tt-0.5)<0.4 && dis4 < disless ) { //only if they are close and the particle is close
			curandState state;

			curand_init(clock64(),1,0,&state);

			double rr1 = curand_uniform(&state);

			if(rr1<probbind) { //bind with rate
			finalbind=mt; //p2 is bound to p1;
			fpos = tt;
		//	cout << p1 << " " << p2 << " " << tt << endl;
			}
		}
}


__global__ void CalculateBindingsGPU(int *list1, int *list2, int *list3, int *list4, float2 *particles, int *bound, double *bound_along, int *changestate,int nt1, int nt2, double L,double l,bool periodic, double disless, double na, double nb, double probbind) {
	//bound is a vector of bindings between particles and the microtubule they are attached too.
	//matrix<double> 


	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index < nt1 ) {
		int p1 = list1[global_index];
		int p2 = list2[global_index];

		int fluid_particle;
		int mt_particle;

		if(p1<p2) {
			fluid_particle = p1;
			mt_particle = p2;
		}
		else{
			fluid_particle = p2;
			mt_particle = p1;
		}
		//printf("global_index: %d, mt_particle: %d, fluid_particle %d\n",global_index,mt_particle,fluid_particle);

		if(bound[fluid_particle]==0&&changestate[fluid_particle]==0) {
			float2 h1 = particles[mt_particle];
			float2 h2 = particles[fluid_particle];
			int init_bind = bound[fluid_particle];
			int fb = bound[fluid_particle];
			double ba = bound_along[fluid_particle];			
		
			float2 uv;
			double dis;
			distance_vector2D(h1,h2,uv,dis, l,periodic);
		// printf("mt_particle %d, fluid_particle %d, (x1,y1) %d,%d (x2,y2) %d,%d\n",mt_particle,fluid_particle,xtemp1,ytemp1,xtemp2,ytemp2);

			BindGPU(mt_particle,fluid_particle,h1,h2,na,nb,L,l,periodic,disless,init_bind,fb,ba, probbind); //only if unbound
			//printf("mt_particle %d, fluid_particle %d, fb %d and ba %f and dis %f\n",mt_particle,fluid_particle,fb,ba, dis);
			if(init_bind!=fb) {
			bound[fluid_particle] =  fb;
			bound_along[fluid_particle] = ba;
			}
		}

			

	}
	else if(global_index >= nt1 && global_index < nt1+nt2) {

		int p1 = list3[global_index-nt1];
		int p2 = list4[global_index-nt1];

		int fluid_particle;
		int mt_particle;

		if(p1<p2) {
			fluid_particle = p1;
			mt_particle = p2;
		}
		else{
			fluid_particle = p2;
			mt_particle = p1;
		}
		//printf("global_index: %d, mt_particle: %d, fluid_particle %d\n",global_index,mt_particle,fluid_particle);

		if(bound[fluid_particle]==0&&changestate[fluid_particle]==0) {
			float2 h1 = particles[mt_particle];
			float2 h2 = particles[fluid_particle];
			int init_bind = bound[fluid_particle];
			int fb = bound[fluid_particle];
			double ba = bound_along[fluid_particle];

			float2 uv;
			double dis;
			distance_vector2D(h1,h2,uv,dis, l,periodic);

			BindGPU(mt_particle,fluid_particle,h1,h2,na,nb,L,l,periodic,disless,init_bind,fb,ba, probbind); //only if unbound
			//printf("mt_particle %d, fluid_particle %d, fb %d and ba %f and dis %f\n",mt_particle,fluid_particle,fb,ba, dis);

			if(init_bind!=fb) {
			bound[fluid_particle] =  fb;
			bound_along[fluid_particle] = ba;
			}
		}	
	}
	else{

	}

	//cout << "calc done" << endl;

}

void Microtubule::callCalculateBindingsGPU(int *list1, int *list2, int *list3, int *list4, float2 *particles, int *bound, double *bound_along, int *changestate, int nt1, int nt2) {
	int blcks = nt1+nt2;
	CalculateBindingsGPU<<<blcks,1>>>(list1,list2,list3,list4,particles,bound,bound_along,changestate,nt1,nt2,Ld,l,is_periodic,disless,(double)na, (double)nb, probbind);
}



template <typename Q>
__global__ void BindingForcesGPUx(float2 *particles, int *bound, double *bound_along, int *con, int *indexf, int *indexp1, int *indexp2, double *forcefx, double *forcefy, double *forcep1x, double *forcep1y,double *forcep2x, double *forcep2y, Q bindp, double L, int na, int nb, double l, bool periodic, int N) {

	//output we want is a list of indices with forces in x ansd y
	int indx1 = threadIdx.x + blockIdx.x * blockDim.x;
	if(indx1> 0 && indx1 < N ) {
		int global_index = con[indx1];
		int mt = bound[global_index];

		if(mt > 0) {
			double tt  =  bound_along[global_index];
			double tt2 = tt*(double)(L-1);
			int k2 = ceil(tt2);
			int k1 = k2 - 1;	

			double t = tt2-k1;

			//na+nb+mt
			int p1a = na+nb+(mt-1)*L+k1;
			int p1b = na+nb+(mt-1)*L+k2;

			//vector1<double> a0(dimension),a1(dimension),a2(dimension);
			float2 a0,a1,a2;

			a0 = particles[global_index];
			a1 = particles[p1a];
			a2 = particles[p1b];

			double dis3=0.0;
			double dis4=0.0;
			float2 uv3,uv4;

			distance_vector2D(a2,a1,uv3,dis3,l,periodic);

			float2 tempvector;
			tempvector.x = a1.x + t*uv3.x;
			tempvector.y = a1.y + t*uv3.y;

			correct_position2D(tempvector,l,periodic);



			distance_vector2D(a0,tempvector,uv4,dis4,l,periodic);


			double f1 = 0;
			
			f1 = bindp(dis4);


			//printf("global index %d and index %d and p1a %d and p1b %d and ba: %f dis3 : %f and dis4 %f and f: %f\n ",indx1,global_index,p1a,p1b,tt,dis3,dis4,f1);

			// for(int j = 0 ; j < dimension ; j++) {
			// forces(p2,j)+=f1*uv4[j]/sqrt(dis4);
			// forces(p1a,j)+=-(1-t)*f1*uv4[j]/sqrt(dis4);
			// forces(p1b,j)+=-(t)*f1*uv4[j]/sqrt(dis4);
			// }
			//int indx1 = con[global_index];
			indexf[indx1-1] = global_index;
			indexp1[indx1-1] = p1a;
			indexp2[indx1-1] = p1b; 

			forcefx[indx1-1] = f1*uv4.x/dis4;

			forcefy[indx1-1] = f1*uv4.y/dis4;

			forcep1x[indx1-1] = -(1-t)*f1*uv4.x/dis4;

			forcep1y[indx1-1] = -(1-t)*f1*uv4.y/dis4;

			forcep2x[indx1-1] = -(t)*f1*uv4.x/dis4;

			forcep2y[indx1-1] = -(t)*f1*uv4.y/dis4;								

		}
	}

}

template<typename Q>
void Microtubule::BindingForcesGPU(float2 *particles, int *bound, double *bound_along, int *&indexf, int *&indexp1, int *&indexp2, double *&forcefx, double *&forcefy, double *&forcep1x, double *&forcep1y,double *&forcep2x, double *&forcep2y, Q bindp, int &nn) {


	int *d_bound;
	int NaNb = na+nb;
	cudaMalloc((void**)&d_bound,(NaNb+1)*sizeof(int));

	setstatefrom<<<NaNb+1,1>>>(d_bound,bound,1,NaNb+1);

	thrust::device_ptr<int> t_bound(d_bound);
	thrust::inclusive_scan(t_bound,t_bound+NaNb+1,t_bound);



	int *d_index;
	cudaMalloc((void**)&d_index,(NaNb+1)*sizeof(int));
	//this->resetindices(d_index,NaNb);
	setstateincr<<<NaNb+1,1>>>(d_index,NaNb+1,-1);

	thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<int> > new_end;

	thrust::device_ptr<int> t_index(d_index);
	new_end = thrust::unique_by_key(t_bound,t_bound+(NaNb+1),t_index);


	//int num_bound = new_end.first;

	int ih1 = thrust::raw_pointer_cast(&new_end.first[0])-thrust::raw_pointer_cast(&t_bound[0])-1;

	cudaMalloc((void**)&indexf,(ih1)*sizeof(int));
	cudaMalloc((void**)&indexp1,(ih1)*sizeof(int));
	cudaMalloc((void**)&indexp2,(ih1)*sizeof(int));

	cudaMalloc((void**)&forcefx,(ih1)*sizeof(double));
	cudaMalloc((void**)&forcefy,(ih1)*sizeof(double));
	cudaMalloc((void**)&forcep1x,(ih1)*sizeof(double));
	cudaMalloc((void**)&forcep1y,(ih1)*sizeof(double));
	cudaMalloc((void**)&forcep2x,(ih1)*sizeof(double));
	cudaMalloc((void**)&forcep2y,(ih1)*sizeof(double));

	// cout << "state set 4" << endl;
	// pausel();
	BindingForcesGPUx<<<ih1+1,1>>>(particles, bound, bound_along, d_index, indexf, indexp1,indexp2,forcefx, forcefy, forcep1x, forcep1y,forcep2x, forcep2y,bindp,  Ld, na, nb, l, is_periodic, ih1+1);


	cudaFree(d_bound);
	cudaFree(d_index);

	nn = ih1;


	// cout << "state set 5" << endl;
	// pausel();	
}

template <typename Q>
void Microtubule::BendingForcesGPU(float2 *particles, int *list1, int *list2, int *list3,double *forcep1x, double *forcep1y, double *forcep2x, double *forcep2y,double *forcep3x, double *forcep3y,Q iny, int &n) {
	calculateforces_threebodyGPU<<<n,1>>>(particles, list1, list2,list3,forcep1x,forcep1y, forcep2x,forcep2y, forcep3x, forcep3y,iny,l, is_periodic, n);
}


__global__ void PositionForcesDueToAnglesGPUx(float2 *particles, int *bound_indices, double *pol, double *v0, int *bound, double *bound_along, int L, double l, bool periodic, double gamma, double dt, int *list1, double *forcep1x, double *forcep1y,int na, int nb, int n) {


	int global_index = threadIdx.x + blockIdx.x * blockDim.x;
	if(global_index > 0 && global_index < n) {
		int indx1 = bound_indices[global_index];

		double tt  =  bound_along[indx1];

		int mt =  bound[indx1];
		double tt2 = tt*(double)(L-1);
		int k2 = ceil(tt2);
		int k1 = k2-1;

		int p1a = na+nb+(mt-1)*L+k1;
		int p1b = na+nb+(mt-1)*L+k2;

		float2 a1,a2;

		a1 = particles[p1a];
		a2 = particles[p1b];

		double dis3;
		float2 uv3;

		distance_vector2D(a1,a2,uv3,dis3,l,periodic);

		uv3.x = uv3.x/dis3;

		uv3.y = uv3.y/dis3;

		double v0_a = v0[global_index];

		double incr = dt*v0_a/gamma;

		double polarity_a = pol[global_index];

		bound_along[indx1] = bound_along[indx1] - polarity_a*incr/((double)(L-1));

	//	printf("global_index: %d and index %d and mt: %d, v0  %d and polarity %d\n",global_index,indx1,mt,v0_a,polarity_a);

		if(bound_along[indx1]<0) bound_along[indx1]=0;

		for(int j = 0 ; j < L ; j++) {
			list1[(global_index-1)*L+j]=na+nb+(mt-1)*L+j;
			forcep1x[(global_index-1)*L+j]=-polarity_a*0.5*v0_a*uv3.x; 
			forcep1y[(global_index-1)*L+j]=-polarity_a*0.5*v0_a*uv3.y;
		}		

		
	}



}

struct polarityfunction {
	int n;
	double pa;
	double pb;

	__host__ __device__ polarityfunction(int nn, double paa, double pbb) : n(nn),pa(paa),pb(pbb) {}
	__host__ __device__ double operator()(int i) {
		if(i<n ) return pa;
		else return pb;
	}
};

void Microtubule::PositionForcesDueToAnglesGPU(float2 *particles, int *bound, double *bound_along, int *&list1, double *&forcep1x, double *&forcep1y, int &n) {
	int *d_bound;
	int NaNb = na+nb;
	cudaMalloc((void**)&d_bound,(NaNb+1)*sizeof(int));

	setstatefrom<<<NaNb+1,1>>>(d_bound,bound,1,NaNb+1);

	thrust::device_ptr<int> t_bound(d_bound);
	thrust::inclusive_scan(t_bound,t_bound+NaNb+1,t_bound);



	int *d_index;
	cudaMalloc((void**)&d_index,(NaNb+1)*sizeof(int));
	//this->resetindices(d_index,NaNb);
	setstateincr<<<NaNb+1,1>>>(d_index,NaNb+1,-1);

	thrust::pair<thrust::device_ptr<int>,thrust::device_ptr<int> > new_end;

	thrust::device_ptr<int> t_index(d_index);
	new_end = thrust::unique_by_key(t_bound,t_bound+(NaNb+1),t_index);


	//int num_bound = new_end.first;

	int ih1 = thrust::raw_pointer_cast(&new_end.first[0])-thrust::raw_pointer_cast(&t_bound[0])-1;


	int ih2 = L*ih1;
	n = ih2;
	if(ih1 > 0) {

	cudaMalloc((void**)&list1,(ih2)*sizeof(int));
	cudaMalloc((void**)&forcep1x,(ih2)*sizeof(double));
	cudaMalloc((void**)&forcep1y,(ih2)*sizeof(double));


	polarityfunction fg(na,polarity_a,polarity_b);
	//setpolarity fg(na,polarity_a,polarity_b);
	// fg.n = na;
	// fg.pa = polarity_a;
	// fg.pb = polarity_b;

	double *d_pol;
	cudaMalloc((void**)&d_pol,(ih1+1)*sizeof(double));
	setstatefunc<<<ih1+1,1>>>(d_pol,d_index,fg,ih1+1);


	polarityfunction fg2(na,v0_a,v0_b);
	// fg2.n = na;
	// fg2.pa = v0_a;
	// fg2.pb = v0_b;	
	
	double *d_v0;
	cudaMalloc((void**)&d_v0,(ih1+1)*sizeof(double));
	setstatefunc<<<ih1+1,1>>>(d_v0,d_index,fg2,ih1+1);





	PositionForcesDueToAnglesGPUx<<<ih1+1,1>>>(particles, d_index, d_pol, d_v0,bound, bound_along, L, l, is_periodic, gamma, dt, list1, forcep1x, forcep1y,na,nb,ih1+1);


	cudaFree(d_v0);
	cudaFree(d_pol);
	}
	cudaFree(d_bound);



}

void prdshft(int &i, int max) {
if(i >= max) i=i-max;
else if(i<0) i=i+max;
else {}
}

void Microtubule::resetchangestate(int *state) {
	int NaNb = na + nb;
	setstate<<<NaNb,1>>>(state,0,na+nb);
}

void Microtubule::resetforce(double *state) {
	double initforce = 0.0;
	setstate<<<totalN,1>>>(state,initforce,totalN);
}


void Microtubule::resetindices(int *index,int n) {
	setstateincr<<<n,1>>>(index,n);
}
