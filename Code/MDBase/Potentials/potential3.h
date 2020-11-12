#ifndef POTENTIAL3_H
#define POTENTIAL3_H

struct potential3 { //this is the base constructor for all potentials.
	bool dl; //does this interaction decay at large distances?
	double interaction_distance;

	

	virtual double energy(vector1<double>&,vector1<double>&)=0; //energy for 3 particles at point 1,2,3
	virtual void force(vector1<double> &, vector1<double> &, vector1<double> &, vector1<double> &, vector1<double> &) = 0;
	//virtual double force2mdx(double)=0; //the factor as a function of distance*squared which when multiplied with dx will give the total force

	virtual void setparameters(const vector1<double>&)=0;
	virtual void printparameters()=0;
	virtual void printparameters(ofstream&)=0;

	virtual potential3* clone() const = 0;
	//virtual potential3* clone() const = 0;
};


struct BendingPotential : potential3 {
	double k0;
	double theta0;

	BendingPotential(double a, double b) : k0(a) , theta0(b) {
		dl = false;
		interaction_distance = 0;
	}

	double energy(vector1<double>&ab,vector1<double>&bc) {
	double a = scalar(ab,bc);
	double mab = sqrt(scalar(ab,ab)); 
	double mbc = sqrt(scalar(bc,bc));

	double theta = acos(a/(mab*mbc));

	return 0.5*k0*SQR(theta-theta0);

	}

	void force(vector1<double> &ab, vector1<double> &bc, vector1<double> &f1, vector1<double> &f2, vector1<double> &f3)
	{ //matrix with force on particle a, particle b and particle c, in order ab is b-a and bc is c-b
		double a = scalar(ab,bc);
		double mabsq = scalar(ab,ab); 
		double mbcsq = scalar(bc,bc);

		double mab = sqrt(mabsq);
		double mbc = sqrt(mbcsq);

		double temp1 = a/(mab*mbc);

		if(temp1>1.) temp1 = 0.999999999;
		else if( temp1<-1 ) temp1 = -0.9999999999;
		else{ }

		double theta = acos(temp1);

		double fac1 = sqrt(1.-SQR(temp1));



		int ds = ab.getsize();

		// matrix<double> forc(3,ab.getsize());

		double fa = -k0*(theta-theta0);
		if(abs(theta-theta0)<1E-10) {
			vector1<double> nullvector(ds);
			f1 = nullvector;
			f2 = nullvector;
			f3 = nullvector;
		}
		else{
			
			
		for(int j = 0 ; j < ds ; j++ )
			f1[j] = fa*(bc[j]*mabsq-ab[j]*a)/(mabsq*fac1*mab*mbc);


		//vector1<double> f2 = (-bc*mabsq*a+(-bc+ab)*mabsq+ab*a*mbcsq)/(mabsq*mbcsq*fac1*sqrt(mabsq)*sqrt(mbcsq));
		for (int j = 0; j < ds; j++)
			f2[j] = fa * ((ab[j] - bc[j]) * mabsq * mbcsq + a * (-(bc[j] * mabsq) + ab[j] * mbcsq)) / (fac1 * mab * mbc * mabsq * mbcsq);

		for (int j = 0; j < ds; j++)
			f3[j] = fa * (bc[j] * a - ab[j] * mbcsq) / (mbcsq * fac1 * mab * mbc );
		}

		// vector1<double> f1 = fa*(bc/(mab*mbc)+ab*a/(SQR(mab)*mab*mbc));
		// vector1<double> f2 = fa*(bc*a/(SQR(mbc)*mbc*mab)-ab*a/(SQR(mab)*mab*mbc)-(bc+ab)/(mab*mbc));
		// vector1<double> f3 = fa*(-bc*a/(SQR(mbc)*mbc*mab)-ab/(mab*mbc));

		// for(int j = 0 ; j < ab.getsize() ; j++ ) {
		// 	forc(0,j)=fa*f1[j];
		// 	forc(1,j)=fa*f2[j];
		// 	forc(2,j)=fa*f3[j];
		// }

		// if(chckmatrixsize(forc,1000)) {
		// 	cout << a << endl;
		// 	cout << mabsq << endl;
		// 	cout << mbcsq << endl;
		// 	cout << fa << endl;
		// 	cout << ab << endl;
		// 	cout << bc << endl;
		// 	cout << f1 << endl;
		// 	cout << f2 << endl;
		// 	cout << f3 << endl;
		// 	error("error in force for bending");
		// }

		//return forc;
	}
	 void setparameters(const vector1<double> &a) {
	 	k0 =  a.gpcons(0);
	 	theta0 =  a.gpcons(1);
	 }
	void printparameters() {
		cout << "harmonic spring constant: " << k0 << " theta0: " << theta0 << endl; 
	}
	void printparameters(ofstream &s) {
		s << k0 << "," << theta0; 
	}	

	BendingPotential* clone() const {
		return new BendingPotential(*this);
	}

};


#endif