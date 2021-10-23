#ifndef NANOTUBE_H
#define NANOTUBE_H

#include "../MDBase/Rotational/LangevinR.h"

//In order to generate the initial sphere conditions, go to "/home/dino/Documents/Nanotube/Sphere.nb"

struct ShellProperties
{
    matrix<int> par;
    matrix<double> posi;
    double k;
    double rm;
};

struct NanotubeAssembly {

int num;

LangevinNVTR *obj;

ComboPatch *pots;

spherical_confinement_3D conf;

// matrix<double> create_spherical_shell(int N, double R);

NanotubeAssembly(double, int); //argument one is radius, argument 2 is the amount

void setviscosity(double);
void setkT(double);

void setpots(ComboPatch &);

void run(int, int, string strbase);

void run_with_real_surface(int,int,ShellProperties &, string strbase);

//AbstractBindingModel *bm;

};

#include "Nanotube.cpp"

#endif /* NANOTUBE_H */
