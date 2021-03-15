#ifndef BIVALENTPATCH_H
#define BIVALENTPATCH_H

struct BivalentPatch : ComboPatch {

    double nx1 = -1.0;
    double ny1 = 0.0;
    double nz1 = 0.0;

    double nx2 = 1.0;
    double ny2 = 0.0;
    double nz2 = 0.0;

    double ang;
    double dis;
    double str;

    // static int *single_patch_ints;
    int *i1;

    BivalentPatch(double strr, double disss, double angg);
    BivalentPatch(const BivalentPatch &patch);

    ~BivalentPatch() { delete i1; }

    int num_patches(const int &);
    void UpdateIterator(const int &i, const int &j);
    void UpdateIteratorSafe(const int &i, const int &j, int **q);
    int get_total_patches(const int &N);
    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj);
    void which_particle(const int &wpi, const int &wpj, int &i, int &j);
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj);

    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);
    
    void CreateFiles();
    
    BivalentPatch *clone() const
    {
        return new BivalentPatch(*this);
    }
};

#include "BivalentPatch.cpp"

#endif /* BIVALENTPATCH_H */
