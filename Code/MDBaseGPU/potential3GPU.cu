struct float2_3 {
	float2 f1,f2,f3;
};


// struct gpupotential3_2D {
// bool dl;
// double interaction_distance;
// __host__ __device__ virtual float2_3 force(float2,float2)=0;

// };

struct BendingPotentialGPU  {
	double k0;
	double theta0;
	bool dl;
	double interaction_distance;

	__host__ __device__ BendingPotentialGPU(double a, double b) : k0(a) , theta0(b), dl(false), interaction_distance(0.0) { }

	__host__ __device__ float2_3 operator()(float2 ab, float2 bc) {
		double a = ab.x*bc.x+ab.y*bc.y;

		double mabsq = ab.x*ab.x+ab.y*ab.y;

		double mbcsq = bc.x*bc.x+bc.y*bc.y;

		double temp1 = a/(sqrt(mabsq*mbcsq));

		if(temp1>1.) temp1 = 0.999999999;
		else if( temp1<-1 ) temp1 = -0.9999999999;
		else{ }

		double theta = acos(temp1);

		double fac1 = sqrt(1.-SQR(temp1));



		double fa = -k0*(theta-theta0);

		float2_3 res;
		if(abs(theta-theta0)<1E-10) {
			
			res.f1.x =0;
			res.f1.y =0;
			res.f2.x =0;
			res.f2.y =0;
			res.f3.x =0;
			res.f3.y =0;
			return res;
		}

		float2 f1;
		float2 f2;
		float2 f3;


		f1.x = (bc.x*mabsq-ab.x*a)/(mabsq*fac1*sqrt(mabsq)*sqrt(mbcsq));
		f1.y = (bc.y*mabsq-ab.y*a)/(mabsq*fac1*sqrt(mabsq)*sqrt(mbcsq));


		//vector1<double> f2 = ((ab - bc)*mabsq*mbcsq + a*(-(bc*mabsq) + ab*mbcsq))/(fac1*sqrt(mabsq*mbcsq)*mabsq*mbcsq);

		f2.x = ((ab.x - bc.x)*mabsq*mbcsq + a*(-(bc.x*mabsq) + ab.x*mbcsq))/(fac1*sqrt(mabsq*mbcsq)*mabsq*mbcsq);
		f2.y = ((ab.y - bc.y)*mabsq*mbcsq + a*(-(bc.y*mabsq) + ab.y*mbcsq))/(fac1*sqrt(mabsq*mbcsq)*mabsq*mbcsq);
		// f2.x = (-bc.x*mabsq*a+(-bc.x+ab.x)*mabsq+ab.x*a*mbcsq)/(mabsq*mbcsq*fac1*sqrt(mabsq)*sqrt(mbcsq));
		// f2.y = (-bc.y*mabsq*a+(-bc.y+ab.y)*mabsq+ab.y*a*mbcsq)/(mabsq*mbcsq*fac1*sqrt(mabsq)*sqrt(mbcsq));

		f3.x = (bc.x*a-ab.x*mbcsq)/(mbcsq*fac1*sqrt(mabsq)*sqrt(mbcsq));
		f3.y = (bc.y*a-ab.y*mbcsq)/(mbcsq*fac1*sqrt(mabsq)*sqrt(mbcsq));

		//printf("%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g,%.20g\n",ab.x,ab.y,bc.x,bc.y,a,mabsq,mbcsq,temp1,theta,fac1,fa,f1.x,f1.y,f2.x,f2.y,f3.x,f3.y);

		res.f1.x = fa*f1.x;
		res.f1.y = fa*f1.y;
		res.f2.x = fa*f2.x;
		res.f2.y = fa*f2.y;
		res.f3.x = fa*f3.x;
		res.f3.y = fa*f3.y;

		return res;
	}
};