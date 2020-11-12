#ifndef COMBOPATCH_H
#define COMBOPATCH_H

#include "potentialtheta.h"




struct ComboPatch {
int **p;
vector1<potentialtheta3D*> potential_bundle; //list of potentials

ComboPatch(int n) : potential_bundle(vector1<potentialtheta3D*>(n)) { }

virtual int num_patches(const int&) = 0; //return the number of patches on particle int
virtual void UpdateIterator(const int&, const int&) = 0; //reassign the pointer to point to something else
virtual int get_total_patches(const int &) = 0;
virtual void which_patch(const int&, const int&, const int &, int&, int & ) = 0; //get patch numbers from particle numbers
virtual void which_particle(const int&, const int&, int&, int&)=0; //get particle numbers from patch numbers
virtual int which_potential(const int&, const int&, const int&, const int& )= 0;
virtual void get_params(const int&,const int&, const int&, int&, int&, int&, int&, int&, int& , double&, double&) = 0;

virtual ComboPatch *clone() const = 0;

};

struct SingPatch : ComboPatch {//single patch
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

    int num_patches(const int&) {return 1;}
    void UpdateIterator(const int &i,const int &j) {
    //do nothing (only one patch)
    }
    int get_total_patches(const int &N) { return N;}
    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj) {
        wpi = i;
        wpj = j;
    }
    void which_particle(const int &wpi, const int &wpj, int &i , int &j ) {
        i = wpi;
        j = wpj;
    }
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj) {
        return 0;
    }
    
    void get_params(const int &i, const int &j, const int &potn, int &nxb1, int &nyb1, int &nzb1, int &nxb2, int &nyb2, int &nzb2, double &d12, double &ang12)
    {
        nxb1  =  nx;
        nyb1  =  ny;
        nzb1  =  nz;
        nxb2  =  nx;
        nyb2  =  ny;
        nzb2  =  nz;
        d12 = dis;
        ang12 =ang;

    }

    SingPatch *clone() const
    {
        return new SingPatch(*this);
    }

};

struct BivalentPatch : ComboPatch {

};


#include "combopatch.cpp"

#endif /* COMBOPATCH_H */
