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
virtual void get_params(const int&,const int&, const int&, double&, double&, double&, double&, double&, double& , double&, double&) = 0;

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

    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12)
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

struct TetrahedralPatch : ComboPatch {

    double nx1 = sqrt(8. / 9.);
    double ny1 = 0.;
    double nz1 = -1. / 3.;

    double nx2 = -sqrt(2. / 9.);
    double ny2 = sqrt(2. / 3.);
    double nz2 = -1. / 3.;

    double nx3 = -sqrt(2. / 9.);
    double ny3 = -sqrt(2. / 3.);
    double nz3 = -1. / 3.;

    double nx4 = 0;
    double ny4 = 0;
    double nz4 = 1.;

    double ang;
    double dis;
    double str;

    int *i1;

    TetrahedralPatch(double strr, double disss, double angg);

    int num_patches(const int &)
    {
        return 4;
    }
    void UpdateIterator(const int &i, const int &j)
    {
        //do nothing (only one patch)
    }
    int get_total_patches(const int &N) { return 4*N; }
    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
    {
        int k1 = potn/4;
        int k2 = potn % 4;
        // i*4 + k1;
        // j*4 + k2;
        wpi = i *4 + k1;
        wpj = j *4 + k2;
    }
    void which_particle(const int &wpi, const int &wpj, int &i, int &j)
    {
        i = wpi/4;
        j = wpj/4;
    }
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
    {
        return (wpi % 4) * 4 + (wpj % 4);   
    }
    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);

    TetrahedralPatch *clone() const
    {
        return new TetrahedralPatch(*this);
    }
};


#include "combopatch.cpp"

#endif /* COMBOPATCH_H */
