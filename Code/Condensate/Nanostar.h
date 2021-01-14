#ifndef NANOSTAR_H
#define NANOSTAR_H

#include "../MDBase/Langevin.h"


struct Nanostar {

    LangevinNVT *obj; //integrator

    double l; //size of the system
    int dimension; //dimension = 3

    int num_nanostars;

    int total_particles;


    int length_of_branch;
    int num_branches;
    //vector<int> stickers;

    vector<mdpair> bindpairs;

    vector<mdtriplet> bendpairs;

    potential *faa;
    potential *hs;
    potential *bindp;
    potential3 *bendp;

    Nanostar(int N, double ll);

    void create_nanostar();


    void Passa_set_nanostar(double theta, double phi, int arms, int armLength, double boxLength, string fileName);

    matrix<double> create_initial_state();
    matrix<double> create_initial_state(string s);

    void set_initial_state(string s);

    matrix<int> gets(matrix<int> &pairs, matrix<int> &specials, matrix<int> &not_specials);

    void run(int total, int every, string);

};

#include "Nanostar.cpp"

#endif /* NANOSTAR_H */
