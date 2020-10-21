#ifndef CONDENSATE_H
#define CONDENSATE_H

#include "../MDBase/LangevinR.h"

struct Condensate {

int num;

LangevinNVTR *obj;

vector1<potentialtheta3D*> potential_bundle;

Condensate(double, int, int);

void set_potential_bundle(vector1<potentialtheta3D*>&);


void run(int,int);



};

#include "Condensate.cpp"
#include "GeometricSchemes.cpp"

#endif
