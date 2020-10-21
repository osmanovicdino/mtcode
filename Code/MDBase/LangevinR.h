#ifndef LANGEVINR_H
#define LANGEVINR_H

#include "Langevin.h"
#include "potentialtheta.h"

class LangevinNVTR : public LangevinNVT {
protected:
	//double gamma;
	vector1<double> im; //MOMENT OF INERTIA TENSOR
	//double dt;
	//double kT;
	//double m;
	//matrix<double> *mom; //momenta of all the particles;
	matrix<double> *angmom; //angular momemntum in 3d
	matrix<double> *orient; //orientation (matrix)
	double gammar;


	double Rt;
	double Rr;
	// double d;//((gamma*dt/2.));
	// double q;//((dt)/2.);
	// double r;//(sqrt(0.5*kT*(gamma)*(m)*(dt)));

	// double c1;// = (dt/m);
	// double c2;// = (1.0/(1.0+(d)));
	// double c3;// = (1.0/(1.0+(d)))*q;
	// double c4;// = (1.0/(1.0+(d)))*r;
	// double c5;// = (1-(d));	
public:
	LangevinNVTR();
	LangevinNVTR(geometry &a); //define a Langevin system with a given geometry

	void initialize(matrix<double> &);
	void initialize(matrix<double> &, matrix<double> &, matrix<double> &, matrix<double> & );

	void setIM(const vector1<double> &);
	void setdt(double);
	void setgamma(double b)  {
		gamma = b;
		Rt = sqrt(2 * gamma * kT / dt);
	}
	void setgammar(double b) 
	{
		gammar = b;
		Rr = sqrt(2 * gammar * kT / dt);
	}

	void setkT(double b) {
		kT = b;
		Rt = sqrt(2 * gamma * kT / dt);
		Rr = sqrt(2 * gammar * kT / dt);
	}

	// genmatx(const vector1<double> &arg);
	// genmaty(const vector1<double> &arg);
	// genmatz(const vector1<double> &arg);
	matrix<double> &getorientation()
	{
		return *(this->orient);
	}
	matrix<double> &getangmom()
	{
		return *(this->angmom);
	}

	double getmag(int i) {
		vector1<double> o = orient->operator()(i);
		return -(o(2) * o(4) * o(6)) + o(1) * o(5) * o(6) + o(2) * o(3) * o(7) -
			o(0) * o(5) * o(7) - o(1) * o(3) * o(8) + o(0) * o(4) * o(8);
	}

	void measured_temperature(ofstream&);

	vector1<double> genfullmat(int);
	void rotate();

	void create_forces_and_torques_sphere(matrix<double> &, matrix<double> &);
	void calculate_forces_and_torques3D(matrix<int> &pairs, potentialtheta3D &, matrix<double> &F, matrix<double> &T);
	void calculate_forces_and_torques3D(matrix<int> &pairs, vector1<potentialtheta3D*> &, matrix<double> &F, matrix<double> &T);
	//matrix<double> calculateforcestheta_pos(matrix<int> &pairs, potentialtheta &);
	//matrix<double> calculateforces_ang(matrix<int> &pairs,potentialtheta&);
	void advancemom_halfstep(matrix<double> &, matrix<double> &);
	void advance_pos();
	// void advancemom_fullstep();
	
	virtual void adv(matrix<int>&);

	//virtual void advm(matrix<int>&,vector1<int>&,vector1<bool>&,intmatrix&); //advance with forces that depend on the internal state of the particle
};

#include "LangevinR.cpp"

#endif