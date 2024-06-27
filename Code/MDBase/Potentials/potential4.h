#ifndef POTENTIAL3_H
#define POTENTIAL3_H

struct potential4 { //this is the base constructor for all potentials.
	bool dl; //does this interaction decay at large distances?
	double interaction_distance;

	virtual double energy(vector1<double> &, vector1<double> &, vector1<double> &, vector1<double> &) = 0; // energy for 3 particles at point 1,2,3
	virtual void force(vector1<double> &, vector1<double> &, vector1<double> &, vector1<double> &) = 0;
	//virtual double force2mdx(double)=0; //the factor as a function of distance*squared which when multiplied with dx will give the total force

	virtual void setparameters(const vector1<double>&)=0;
	virtual void printparameters()=0;
	virtual void printparameters(ofstream&)=0;

	virtual potential4* clone() const = 0;
	//virtual potential3* clone() const = 0;
};


};


#endif