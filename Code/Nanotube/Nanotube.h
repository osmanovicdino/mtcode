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

    void DoAnMC(double, bool);
};

struct particle_adder {
int No_types;
vector1<int> nums;
vector1<int> pats;



};

struct NanotubeAssembly {

int num;

LangevinNVTR *obj;

TetrahedralWithBivalent *pots;

spherical_confinement_3D conf;

double myrmax;

double ll;

// matrix<double> create_spherical_shell(int N, double R);

NanotubeAssembly(double, int); //argument one is radius, argument 2 is the amount

void setviscosity(double);
void setkT(double);

void add_particle2();

void add_particle42(int);

void setpots(TetrahedralWithBivalent &);

void run(int, int, string strbase);

void run_with_real_surface(int, int, ShellProperties &, string strbase);
void run_with_real_surface_add_particles(int, int, ShellProperties &, double prod, string strbase);

//AbstractBindingModel *bm;

};

#include "Nanotube.cpp"

#endif /* NANOTUBE_H */
