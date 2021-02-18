#ifndef CONDENSATE_H
#define CONDENSATE_H

#include "../MDBase/Rotational/LangevinR.h"

struct Condensate {

int num;

LangevinNVTR *obj;

ComboPatch *pots;

AbstractBindingModel *bm;

Condensate(double, int);

void setviscosity(double);
void setkT(double);

void setpots(ComboPatch&);

void setBindingModel(AbstractBindingModel&);

void run(int,int, string strbase);

void run_singlebond(int,int,string strbase);

void run_singlebond_different_sizes(int runtime, int every, int div, string strbase);

void run_singlebond_different_sizes_continue(int runtime, int every, int div, int startval, BinaryBindStore &bbs2, string strbase);
};

#include "Condensate.cpp"
//#include "GeometricSchemes.cpp"

#endif
