#ifndef SINGPATCH_H
#define SINGPATCH_H

struct SingPatch : ComboPatch
{                                 //single patch
    int patches_per_particle = 1; //as the number of patches don't matter;

    double nx = 1.0;
    double ny = 0.0;
    double nz = 0.0;

    double ang;
    double dis;
    double str;

    // static int *single_patch_ints;
    int *i1;

    SingPatch(double strr, double disss, double angg);

    ~SingPatch() { delete i1; }

    int num_patches(const int &);
    void UpdateIterator(const int &i, const int &j);
    void UpdateIteratorSafe(const int &i, const int &j, int **q);
    int get_total_patches(const int &N);
    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj);
    void which_particle(const int &wpi, const int &wpj, int &i, int &j);
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj);

    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);
    SingPatch *clone() const
    {
        return new SingPatch(*this);
    }
};

#include "SingPatch.cpp"

#endif /* SINGPATCH_H */
