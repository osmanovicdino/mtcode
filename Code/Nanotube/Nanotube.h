#ifndef NANOTUBE_H
#define NANOTUBE_H

#include "../MDBase/Rotational/LangevinR.h"
#include "ParticleAdder.h"

struct WeiM
{
    double weight;
    int M;
};
//In order to generate the initial sphere conditions, go to "/home/dino/Documents/Nanotube/Sphere.nb"

struct ShellProperties
{
    matrix<int> par;
    matrix<double> posi;
    double k;
    double rm;

    void DoAnMC(double, bool);
};


struct NanotubeAssembly {

int num;

LangevinNVTR *obj;

ComboPatch *pots;

spherical_confinement_3D conf;

double myrmax;

double ll;

// matrix<double> create_spherical_shell(int N, double R);

NanotubeAssembly(double, int); //argument one is radius, argument 2 is the amount

void setviscosity(double);
void setkT(double);

void add_particle2();

void add_particle42(int);

void setpots(ComboPatch &);

matrix<double> calculate_covariance(int Ns);

void run(int, int, string strbase);

void run_anneal(int, int, int, string strbase);

void run_with_real_surface(int, int, ShellProperties &, matrix<double> &constantF, string strbase);
void run_with_real_surface_add_particles(int, int, ShellProperties &, double prod, WeiM &c1, string strbase);
void run_with_real_surface_add_particles_continue(int, int, int, ShellProperties &, double prod, WeiM &c1, matrix<double>&, matrix<double>&, vector1<int>&,string strbase);
void run_add_particles(int, int, double prod, string strbase);

//AbstractBindingModel *bm;

};

#include "Nanotube.cpp"

#endif /* NANOTUBE_H */
