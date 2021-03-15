#ifndef NANOTUBE_H
#define NANOTUBE_H

#include "../MDBase/Rotational/LangevinR.h"

struct NanotubeAssembly {

int num;

LangevinNVTR *obj;

ComboPatch *pots;

spherical_confinement_3D conf;

NanotubeAssembly(double, int); //argument one is radius, argument 2 is the amount

void setviscosity(double);
void setkT(double);

void setpots(ComboPatch &);

void run(int, int, string strbase);

//AbstractBindingModel *bm;

};

#include "Nanotube.cpp"

#endif /* NANOTUBE_H */
