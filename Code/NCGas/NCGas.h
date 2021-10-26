#ifndef NCGAS_H
#define NCGAS_H



struct NCGas {
LangevinNVT *obj;
int num;
double eps;
double eqeps;
double lambda;
double l;

double sigma;
double eta;

NCGas(double, int, double, double );


double chm;
double lam2;
double lam4;

void seteps(double);
void setkT(double);
void seteqeps(double, double);
void setgeo(cube&);

matrix<double> calculate_virial_with_matrices(matrix<double>&,matrix<double>&);
matrix<double> calculate_virial();

void run(int);




//vector1<double> forcecalc1(double,double,double,double,double,double);
void forcecalc2(int, int, vector1<double>&);
//vector1<double> forcecalc3(int);

};

#include "NCGas.cpp"

#endif