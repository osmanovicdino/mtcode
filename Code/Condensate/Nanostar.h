#ifndef NANOSTAR_H
#define NANOSTAR_H

#include "../MDBase/Langevin.h"


struct Nanostar {

    LangevinNVT *obj; //integrator

    double l; //size of the system
    int dimension; //dimension = 3

    int num_nanostars;

    int total_particles;


    
    int num_branches;
    vector1<int> length_of_each_branch;

    vector1<int> accumulate;

    int total_parts_per_nanostar;
    //vector<int> stickers;

    vector<mdpair> bindpairs;

    vector<mdtriplet> bendpairs;

    vector<int> stickerList;

    potential *faa;
    potential *hs;
    potential *bindp;
    potential3 *bendp;


    Nanostar(int N, double ll, vector1<int>&, matrix<double> &);

    void create_nanostar();


    
    matrix<double> create_initial_state(matrix<double> &ori);
    matrix<double> create_initial_state(string s);

    void set_initial_state(string s);

    void gets(matrix<int> &pairs, matrix<int> &specials, matrix<int> &not_specials);

    void DoAnMC();

    void run(int total, int every,int, string);
    
};

#include "Nanostar.cpp"

#endif /* NANOSTAR_H */
