#ifndef NANOTUBE_H
#define NANOTUBE_H

#include "../MDBase/Rotational/LangevinR.h"
#include "ParticleAdder.h"
#include "GeneticDataStructure.h"

struct WeiM
{
    double weight;
    int M;
};
//In order to generate the initial sphere conditions, go to "/home/dino/Documents/Nanotube/Sphere.nb"

struct ShellProperties
{
    matrix<int> par;
    matrix<int> quad;
    matrix<double> posi;
    matrix<double> bindingdis;
    double k;
    double rm;
    double kappa;

    ShellProperties(matrix<int> parr,matrix<int> quadd,matrix<double> posii, matrix<double> bindingdiss,double kk, double rmm, double kappaa) : par(parr), quad(quadd), posi(posii),bindingdis(bindingdiss), k(kk), kappa(kappaa), rm(rmm) {

    }

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
NanotubeAssembly(double,int,bool); //just add randomly

void setviscosity(double);
void setkT(double);

void add_particle2();

void add_particle42(int);

void setpots(ComboPatch &);

matrix<double> calculate_covariance(int Ns);

void run(int, int, string strbase);

void run_anneal(int, int, int, string strbase);

void run_with_real_surface(int, int, ShellProperties &, matrix<double> &constantF, string strbase, bool);
void run_with_real_surface_add_particles(int, int, ShellProperties &, double prod, WeiM &c1, string strbase);
void run_with_real_surface_add_particles_continue(int, int, int, ShellProperties &, double prod, WeiM &c1, matrix<double>&, matrix<double>&, vector1<int>&,string strbase);
void run_with_real_surface_remove_add_particles_continue(int, int, int, ShellProperties &, double prod, double, double, WeiM &c1, matrix<double> &, matrix<double> &, vector1<int> &, string strbase);
void run_add_particles(int, int, double prod, string strbase);
void run_bending_modulus(int, int, double,string);

void run_with_real_surface_add_particles_steric(int, int, ShellProperties &, double prod, WeiM &c1, string strbase);

void run_no_patchy(int, int, ShellProperties &, double prod, double ks, double ks2, string strbase);

void run_no_patchy_noshell(int, int, double prod, double ks, double ks2, string strbase);

void run_no_patchy_continue(int, int, int, ShellProperties &, double prod, double ks, double ks2, matrix<double>&, vector1<int>&, vector<vector<int>>&, string strbase);

void run_box(int,int,double,geneticcode&,string);
void run_box_equil(int, int, double, geneticcode &, string);
void run_box_equil_cont(int, int, double, geneticcode &, string, matrix<double>&,matrix<double>&, matrix<int>&);
void run_box_anchors(int, int, double, int, geneticcode &, string);
//AbstractBindingModel *bm;

};

#include "Nanotube.cpp"
#include "Nanotube2.cpp"
#include "Nanotube3.cpp"

#endif /* NANOTUBE_H */
