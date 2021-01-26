#ifndef NANOSTAR_H
#define NANOSTAR_H

#include "../MDBase/Langevin.h"
#include <vector>
#include "../MathHelper/Helper.h"

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

    vector<mdtriplet> bendtriples;
    matrix<double> particles;
    vector<int> stickerList;

    potential *faa;
    potential *hs;
    potential *bindp;
    potential3 *bendp;

    Nanostar(int N, double ll);

    void create_nanostar();
    vector1<double> initCoord;

    void Passa_set_nanostar(vector1<double> start, double theta, double phi, int arms, int armLength, string fileName);
    void setNanostar(vector1<double> start, int arms, int armLength, double theta, string fileName);
    void sortPairsTriplets(int arms, int armLength);
    void initStickerList(int arms, int armLength);
    void inStickerList(matrix<int> &possiblePairs, matrix <int> &specials, matrix<int> &notSpecials);
    matrix<double> create_initial_state();
    matrix<double> create_initial_state(string s);

    void set_initial_state(string s);

    matrix<int> gets(matrix<int> &pairs, matrix<int> &specials, matrix<int> &not_specials);

    void run(int total, int every, string);

};

#include "Nanostar.cpp"

#endif /* NANOSTAR_H */
