#ifndef NCGASR_H
#define NCGASR_H

#include "Langevin.h"


struct NCGasR {
LangevinNVT *obj;
int num;
double eps;
double eqeps;
double lambda;
double l;
matrix<double> *orient;

double angd;
double xi;
double v0;

double fac;
double dt;

// double exttheta;
// double extphi;
// double ext_strength;

NCGasR(double, int);


double chm;
double lam2;
double lam4;

void seteps(double);
void setkT(double);
void seteqeps(double);
void setgeo(cube&);

matrix<double> PositionForcesDueToAngles_mips(); //for example mips

matrix<double> AngleForcesDueToPosition(matrix<int> &pairs);

matrix<double> AngleForcesOnly(double,double,double);

void UpdateAngles(matrix<double>&);

void run(int);




//vector1<double> forcecalc1(double,double,double,double,double,double);
//void forcecalc2(int, int, vector1<double>&);
//vector1<double> forcecalc3(int);

};

#include "NCGasR.cpp"

#endif