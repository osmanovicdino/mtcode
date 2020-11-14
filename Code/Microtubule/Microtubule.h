
#ifndef MICROTUBULE_H
#define MICROTUBULE_H

struct float2;
struct gpupotential;
struct gpupotential3_2D;

struct Microtubule {
LangevinNVT *obj;


int num;
int dimension;
// double eps;
// double eqeps;
// double lambda;
double l;
bool is_periodic;
// double v00;

int na; //number of type A
int nb; //number of type B
int nc; //number of type C
int number_of_microtubules; // number of microtubules
int totalN; //totalN
matrix<double> *orient; // orientation of microtubules
int L; //Length of the microtubles
double Ld;
double probunbind;
double probbind;

bool spatial_dependence;
double gamma;
double max_s; //the maximum speed the microtubule can go;


double v0_a;
double v0_b;

double disless;

double polarity_a;
double polarity_b;

double excess_force_distance;



vector1<int> pai; //indices of particles of type A
vector1<int> pbi; //indices of particles of type B
vector1<int> pci; //indices of particles of type C

//matrix<bool> binding; //matrix of bindings
vector1<int> bound; //is it bound for fluid
vector1<double> bound_along;
//vector1<int> micrtbound; //is it bound microtubule

potential *faa;
potential *fbb;
potential *fcc;
potential *fab;
potential *fbc;
potential *fac;
potential *bindp;
potential *bindm;
potential3 *bendp;

matrix<int> *bondpairs;
matrix<int> *bendtriplets;


double dt;

// double exttheta;
// double extphi;
// double ext_strength;

Microtubule(double, int, int, int, int);
~Microtubule();

//void seteps(double);
void setkT(double);
//void seteqeps(double);
void setgeo(cube&);
void setgamma(double);

void setpotaa(potential&);
void setpotab(potential&);
void setpotac(potential&);
void setpotbb(potential&);
void setpotbc(potential&);
void setpotcc(potential&);
// void setv00(double);

void setprobunbind(double);
void setprobbind(double);
void set_excess_force_distance(double);
void setv0(double, double);
void setpolarity(double, double);

void set_initial_conditions(string filename1, string filename2, string filename3);

void set_gamma_spatial_dependence(bool);

void resetchangestate(int *);

void resetindices(int *, int);

void resetforce(double *);

matrix<double> TubeForces(int&,int&,int&);
void Bind(int&,int&,int,int&,double&);

// void CalculateUnbindingsGPU();
// void CalculateBindingsGPU();
//void callBindGPU(int &, int &, float2, float2, double , double , bool , double, int , int&, double &);

void callCalculateUnbindingsGPU(float2 *, int *, double *, int *);

void callCalculateBindingsGPU(int *, int *, int *, int *, float2 *, int*, double*, int *,int , int );

void CalculateBindings(matrix<int>&,matrix<int>&);

matrix<double> BindingForces();

template <typename Q>
void BindingForcesGPU(float2 *, int *, double *, int *&, int *&, int *&, double *&, double *&, double *&, double *&,double *&, double *&, Q, int &);

template <typename Q>
void BendingForcesGPU(float2 *, int *, int *, int *,double *, double *, double *, double *,double *, double *, Q, int &);

matrix<double> BindingForcesTest(vector1<int>&,vector1<double>&,matrix<double>&);

matrix<double> constantMTforce(); //propulsion force

matrix<double> constantMTforce(vector1<int>&,vector1<int>&); //propulsion force, in the direction of particles i1 and i2

matrix<double> PositionForcesDueToAngles(); //propulsion force

void PositionForcesDueToAnglesGPU(float2 *, int *, double *, int *&, double *&, double *&, int&);


void minimal_distance_between_point_and_line(vector1<double> &point1, vector1<double> &end1, vector1<double> &end2, double&,vector1<double>&,double&);

void minimal_distance_between_point_and_line(int&, int&, int&, double&,vector1<double>&,double&);
//matrix<double> PositionForcesDueToAngles(double); //for example mips

void ForcesDueToPositionPL(matrix<int> &pairs,matrix<double>&); //interaction of angular forces, i.e., if the force is not between the centres of the rods
void ForcesDueToPositionLL(matrix<int> &pairs,matrix<double>&,matrix<double>&); //interaction of angular forces, i.e., if the force is not between the centres of the rods


void BendingForces(matrix<int>); //matrix of triplets
//matrix<double> AngleForcesOnly(double,double,double);

void UpdateAngles(matrix<double>&,matrix<double>&);

void UpdateBoundAlong();

// void run2D(int);

void run(int,int);


template <typename Fun>
void runSV(int,int,Fun);


template <typename Fun>
void runPV(int,int,Fun);

void runGPU(int,int);

template <typename Fun>
void runGPUSV(int,int,Fun);

template <typename Fun>
void runGPUPV(int,int,Fun);

template <typename Fun>
void runGPUcheck(int,int,Fun);

template <typename Fun>
void runMTONLY(int,int,Fun);

template <typename Fun>
void runMTONLY_initialstate(int, int, Fun, double, vector1<int>, vector1<int>);

//vector1<double> forcecalc1(double,double,double,double,double,double);
//void forcecalc2(int, int, vector1<double>&);
//vector1<double> forcecalc3(int);

};

#include "Microtubule.cpp"
#include "MicrotubuleRuns.cpp"

// #include "../MicrotubuleGPU/MicrotubuleGPU.cu"

// #include "../MicrotubuleGPU/MicrotubuleGPUrun.cu"

#endif