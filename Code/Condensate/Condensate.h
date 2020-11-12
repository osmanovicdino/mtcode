#ifndef CONDENSATE_H
#define CONDENSATE_H

#include "../MDBase/Rotational/LangevinR.h"

struct Condensate {

int num;

LangevinNVTR *obj;

ComboPatch *pots;

Condensate(double, int);

void setpots(ComboPatch*);

void run(int,int, string strbase);



};

#include "Condensate.cpp"
//#include "GeometricSchemes.cpp"

#endif
