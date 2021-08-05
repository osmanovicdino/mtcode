#ifndef CONDENSATE_H
#define CONDENSATE_H

#include "../MDBase/Rotational/LangevinR.h"

struct Condensate {

int num;
double ls;
double size_mol = 1.0;

LangevinNVTR *obj;

ComboPatch *pots;

AbstractBindingModel *bm;

Condensate(double, int);

void setviscosity(double);
void setkT(double);

void setpots(ComboPatch&);

void setBindingModel(AbstractBindingModel&);

void setup_tight_packing(double size);

void setup_large_droplet(int N1, int N2, int N3, double dens, double ll);



void run(int,int, string strbase);

void run_singlebond(int,int,string strbase);

void run_singlebond_continue(int, int, int startval, BinaryBindStore &bbs2, string strbase);

void run_singlebond_eq(int,int,string strbase);

void run_singlebond_different_sizes(int runtime, int every, int div, double size1, double size2, string strbase);

void run_singlebond_different_sizes_continue(int runtime, int every, int div, double size1, double size2, int startval, BinaryBindStore &bbs2, string strbase);

void run_singlebond_different_sizes_continue_thetalist(int runtime, int every, int div, double size1, double size2, int startval, BinaryBindStore &bbs2, string strbase);


};

#include "Condensate.cpp"
//#include "GeometricSchemes.cpp"

#endif
