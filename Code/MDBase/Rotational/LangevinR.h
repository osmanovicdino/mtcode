#ifndef LANGEVINR_H
#define LANGEVINR_H

#include "../Langevin.h"
#include "../Potentials/potentialtheta.h"
#include "../Bindings/GraphAlgorithms.h"
#include "../Bindings/BinaryBindStore.h"
#include "../Bindings/BindingModelFull.h"
#include "../Bindings/BindingModelSingle.h"
#include "../Bindings/BindingModelBinary.h"
#include "../Bindings/BindingModelTernary.h"
#include "../Potentials/ComboPatch/combopatch.h"

struct thetapair
{
	double thetai;
	double thetaj;
};

struct patchint { //for storing information about all possible pairs
	unsigned int particle_index1;
	unsigned int particle_index2;
	unsigned int patch_index1;
	unsigned int patch_index2;
	unsigned int potn;

	inline void get_all(int &a1, int &a2, int &p, int &pi1, int &pi2) {
		a1= particle_index1;
		a2= particle_index2;
		p=potn;
		pi1=patch_index1;
		pi2=patch_index2;
	}

	inline void get_particle(int &a1, int &a2) {
		a1 = particle_index1;
		a2 = particle_index2;
	}

	inline void get_patch(int &p, int &pi1, int &pi2) {
		p = potn;
		pi1 = patch_index1;
		pi2 = patch_index2;
	}

	friend inline bool check_same_particle(const patchint &a1, const patchint &a2) {
		if(a1.particle_index2 == a2.particle_index2) {
			if(a1.particle_index1 == a2.particle_index1) {
				return true;
			}
			return false;
		}
		else {
			return false;
		}
	}
};

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
	LangevinNVTR(cube &a); //define a Langevin system with a given geometry
	LangevinNVTR(const LangevinNVTR&);
	LangevinNVTR& operator=(const LangevinNVTR&);
	LangevinNVTR *clone() const { 
		return new LangevinNVTR(*this); }
	~LangevinNVTR();

	void initialize(matrix<double> &);
	void initialize(matrix<double> &, matrix<double> &, matrix<double> &, matrix<double> & );

	void printall() {
		cout << "printing all parameters" << endl;
		cout << angmom->getNsafe() << endl;
		cout << orient->getNsafe() << endl;
		cout << mom->getNsafe() << endl;
		cout << dat->getNsafe() << endl;
		cout << "printed" << endl;
		cout << endl;

	}

	void setIM(const vector1<double> &);
	void setdt(double);
	void setgamma(double b)  {
		gamma = b;
		d = (gamma * dt / 2.);
		q = (dt) / 2.;
		r = sqrt(0.5 * kT * (gamma) * (m) * (dt));

		c1 = (dt / m);
		c2 = (1.0 / (1.0 + (d)));
		c3 = (1.0 / (1.0 + (d))) * q;
		c4 = (1.0 / (1.0 + (d))) * r;
		c5 = (1 - (d));

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
	void setorientation(matrix<double> &ort)  {
		// cout << ort.getNsafe() << endl;
		// matrix<double> *res = ort.clone();
		delete orient;
		orient = ort.clone();
		//orient = res;

		int jdimension = orient->getncols();
		if (jdimension != 3*(geo.dimension) )
		{
			cout << dimension << " " << geo.dimension << endl;
			error("set orient dimensions must match in MD");
		}
	}

	void setangularmomenta(matrix<double> &amom) {
		//if(amom.getNsafe() != angmom->getNsafe() || amom.getncols() != angmom->getncols() ) error("not setting correct dimensions for setting angular momentum");
		
		delete angmom;
		

		angmom =  amom.clone();
	}

	void setup_random_binding(vector<patchint> &pairs, vector<int> &divs, ComboPatch &iny, BinaryBindStore &bo, AbstractBindingModel &bm);


	void setorientation(vector1<double> &a, vector1<double> &b, int i) {
		//maps a to b
		// (*orient)(i,0) = 1 - (SQR(a(1)*b(0) - a(0)*b(1) + SQR(a(2)*b(0) - a(0)*b(2))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		// (*orient)(i,1) = -(a(1)*b(0)) + a(0)*b(1) + ((-(a(2)*b(0)) + a(0)*b(2))*(a(2)*b(1) - a(1)*b(2)))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		// (*orient)(i,2) = -(a(2)*b(0)) + a(0)*b(2) + ((-(a(1)*b(0)) + a(0)*b(1))*(-(a(2)*b(1)) + a(1)*b(2)))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		// (*orient)(i,3) = a(1)*b(0) - a(0)*b(1) + ((a(2)*b(0) - a(0)*b(2))*(-(a(2)*b(1)) + a(1)*b(2)))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		// (*orient)(i,4) = 1 - ((SQR(a(0) + SQR(a(2))*SQR(b(1) - 2*a(1)*b(1)*(a(0)*b(0) + a(2)*b(2)) + SQR(a(1)*(SQR(b(0) + SQR(b(2)))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		// (*orient)(i,5) = -(a(2)*b(1)) + a(1)*b(2) + ((a(1)*b(0) - a(0)*b(1))*(-(a(2)*b(0)) + a(0)*b(2)))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		// (*orient)(i,6) = a(2)*b(0) - a(0)*b(2) + ((a(1)*b(0) - a(0)*b(1))*(a(2)*b(1) - a(1)*b(2)))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		// (*orient)(i,7) = a(2)*b(1) - a(1)*b(2) + ((-(a(1)*b(0)) + a(0)*b(1))*(a(2)*b(0) - a(0)*b(2)))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		// (*orient)(i,8) = 1 - (SQR(a(2)*(SQR(b(0) + SQR(b(1)) - 2*a(2)*(a(0)*b(0) + a(1)*b(1))*b(2) + (SQR(a(0) + SQR(a(1))*SQR(b(2))/(1 + a(0)*b(0) + a(1)*b(1) + a(2)*b(2));
		(*orient)(i, 0) = 1 - (SQR(a[1] * b[0] - a[0] * b[1]) + SQR(a[2] * b[0] - a[0] * b[2])) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		(*orient)(i, 1) = -(a[1] * b[0]) + a[0] * b[1] + ((-(a[2] * b[0]) + a[0] * b[2]) * (a[2] * b[1] - a[1] * b[2])) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		(*orient)(i, 2) = -(a[2] * b[0]) + a[0] * b[2] + ((-(a[1] * b[0]) + a[0] * b[1]) * (-(a[2] * b[1]) + a[1] * b[2])) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		(*orient)(i, 3) = a[1] * b[0] - a[0] * b[1] + ((a[2] * b[0] - a[0] * b[2]) * (-(a[2] * b[1]) + a[1] * b[2])) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		(*orient)(i, 4) = 1 - ((SQR(a[0]) + SQR(a[2])) * SQR(b[1]) - 2 * a[1] * b[1] * (a[0] * b[0] + a[2] * b[2]) + SQR(a[1]) * (SQR(b[0]) + SQR(b[2]))) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		(*orient)(i, 5) = -(a[2] * b[1]) + a[1] * b[2] + ((a[1] * b[0] - a[0] * b[1]) * (-(a[2] * b[0]) + a[0] * b[2])) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		(*orient)(i, 6) = a[2] * b[0] - a[0] * b[2] + ((a[1] * b[0] - a[0] * b[1]) * (a[2] * b[1] - a[1] * b[2])) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		(*orient)(i, 7) = a[2] * b[1] - a[1] * b[2] + ((-(a[1] * b[0]) + a[0] * b[1]) * (a[2] * b[0] - a[0] * b[2])) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		(*orient)(i, 8) = 1 - (SQR(a[2]) * (SQR(b[0]) + SQR(b[1])) - 2 * a[2] * (a[0] * b[0] + a[1] * b[1]) * b[2] + (SQR(a[0]) + SQR(a[1])) * SQR(b[2])) / (1 + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
	}
	
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

	vector1<double> transformvector(vector1<double> &v, int p1) {
		double qtemp0 = orient->gpcons(p1, 0);
		double qtemp1 = orient->gpcons(p1, 1);
		double qtemp2 = orient->gpcons(p1, 2);
		double qtemp3 = orient->gpcons(p1, 3);
		double qtemp4 = orient->gpcons(p1, 4);
		double qtemp5 = orient->gpcons(p1, 5);
		double qtemp6 = orient->gpcons(p1, 6);
		double qtemp7 = orient->gpcons(p1, 7);
		double qtemp8 = orient->gpcons(p1, 8);

		vector1<double> newv(3);

		double nxb1 = v[0];
		double nyb1 = v[1];
		double nzb1 = v[2];

		newv[0] = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
		newv[1] = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
		newv[2] = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;
		return newv;
	}
	void measured_temperature();
	
	void measured_temperature(ofstream&);
	void measured_temperature(vector1<double>&);

	vector1<double> genfullmat(int);
	void rotate();

	template <class vec>
	void rotate(vec&);

	// void one_bond_per_patch_condition(matrix<int> &pairs, vector<potentialtheta3D*> &);
	//void create_random_forces(matrix<double>&, matrix<double>&);
	template <class vec>
	void create_forces_and_torques_sphere(matrix<double> &, matrix<double> &, matrix<double> &, vec&, bool);
	
	void create_forces_and_torques_sphere(matrix<double> &, matrix<double> &, matrix<double> &); //transform between body fixed frame and lab fixed frame
	void calculate_forces_and_torques3D(matrix<int> &pairs, potentialtheta3D &, matrix<double> &F, matrix<double> &T); //calculation of force for a single angle dependent theta
	void calculate_forces_and_torques3D(matrix<int> &pairs, vector1<potentialtheta3D*> &, matrix<double> &F, matrix<double> &T); //calculation of force for a potential bundle
	void calculate_forces_and_torques3D(matrix<int> &pairs, ComboPatch &, matrix<double> &F, matrix<double> &T); //calculation of force for a potential bundle
	void calculate_forces_and_torques3D(vector<patchint> &pairs, ComboPatch &, matrix<double> &F, matrix<double> &T); // calculation of force for a potential bundle

	void calculate_forces_and_torques3D_onlyone(matrix<int> &pairs, vector1<potentialtheta3D *> &, BinaryBindStore & , AbstractBindingModel&, matrix<double> &F, matrix<double> &T); //calculation of force for a potential bundle
	matrix<int> CreateEdgeList(matrix<int> &adj, vector1<int> &lens);
	void calculate_forces_and_torques3D_onlyone(matrix<int> &pairs, ComboPatch&, BinaryBindStore &, AbstractBindingModel &, matrix<double> &F, matrix<double> &T); //calculation of force for a potential bundle
	void calculate_forces_and_torques3D_onlyone_nonlets(matrix<int> &pairs, ComboPatch &, BinaryBindStore &, AbstractBindingModel &, matrix<double> &F, matrix<double> &T); //calculation of force for a potential bundle
	void calculate_forces_and_torques3D_onlyone_nonlets(vector<patchint> &pairs, vector<int> &divs, ComboPatch &, BinaryBindStore &, AbstractBindingModel &, matrix<double> &F, matrix<double> &T); //calculation of force for a potential bundle
	void calculate_forces_and_torques3D_onlyone_nonlets_eq(vector<patchint> &pairs, vector<int> &divs, ComboPatch &, BinaryBindStore &, AbstractBindingModel &, matrix<double> &F, matrix<double> &T); //calculation of force for a potential bundle

	//matrix<double> calculateforcestheta_pos(matrix<int> &pairs, potentialtheta &);
	//matrix<double> calculateforces_ang(matrix<int> &pairs,potentialtheta&);
	void advancemom_halfstep(matrix<double> &, matrix<double> &);

	template <class vec>
	void advancemom_halfstep(matrix<double> &, matrix<double> &, vec&);
	//void advance_pos();
	// void advancemom_fullstep();
	vector<thetapair> check_arg_thetas_per_pair(matrix<int> &pairs, ComboPatch &iny);
	vector<patchint> calculate_patch_list(matrix<int> &pairs, ComboPatch &iny);
	vector<patchint> calculate_bound_pairs(matrix<int> &pairs, ComboPatch &iny);
	virtual void adv(matrix<int>&);

	//virtual void advm(matrix<int>&,vector1<int>&,vector1<bool>&,intmatrix&); //advance with forces that depend on the internal state of the particle
};

#include "LangevinR.cpp"
#include "LangevinRforceUtil.cpp"
#include "LangevinRforce.cpp"

#endif