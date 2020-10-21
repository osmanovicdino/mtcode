// struct gpupotential {
// bool dl;
// double interaction_distance;
// __host__ __device__ double force(double x) {
// 	printf("error, base constructor method called%d\n",(int)x);
// 	}
// };
/* WRITTEN WITHOUT VIRTUAL BASE CLASS AS CUDA DOES NOT SUPPORT IT*/


struct WCApotentialGPU {
double epsilon; //epsilon
double sigma; //sigma
double attrepsilon; //epsilon modifier
double paramlimit;
double sigma6;
double sigma12;
bool dl;
double interaction_distance;

	__host__ __device__ WCApotentialGPU(double a, double b, double c) : epsilon(a),sigma(b),attrepsilon(c),sigma6(b*b*b*b*b*b),sigma12(b*b*b*b*b*b*b*b*b*b*b*b),interaction_distance(2.5*b),dl(true),paramlimit(1.12246*b) {
		// sigma6= SQR(CUB(sigma)); sigma12 = SQR(sigma6);
		// interaction_distance=2.5*sigma;
		// dl = true;
		// paramlimit=1.12246*sigma;
	}


	__host__ __device__ double operator()(double x) {
	    if(x < paramlimit) {
	    double a;
	    a=x*x*x;
	    
	    a=SQR(a);

		return 4*epsilon*(12*sigma12/(x*SQR(a))-6*sigma6/(x*a));
		}
		else if(x >= paramlimit && x< interaction_distance) {
	    double a;
	    a=x*x*x;
	    
	    a=SQR(a);

		return 4*attrepsilon*(12*sigma12/(x*SQR(a))-6*sigma6/(x*a));			
		}
		else {
			return 0;
		}
	}

};

struct HarmonicPotentialGPU {
	double k;
	double x0;
	bool dl;
	double interaction_distance;	

	__host__ __device__ HarmonicPotentialGPU(double a, double b) : k(a),x0(b),dl(false),interaction_distance(0.0) { }


	__host__ __device__ double operator()(double x) {
		return -k*(x-x0);
	}

};

struct FENEPotentialGPU {
	double kbond;
	double R_0;
	bool dl;
	double interaction_distance;	

	__host__ __device__ FENEPotentialGPU(double a, double b) : kbond(a),R_0(b),dl(false),interaction_distance(0.0) { }


	__host__ __device__ double operator()(double x) { return -kbond*x/(1-SQR(x)/SQR(R_0));}


};