#ifndef COMBOPATCH_H
#define COMBOPATCH_H

#include "potentialtheta.h"




struct ComboPatch {
int **p; //pointer to array where potentials are
bool safe; // is this combo patch parallel safe?
vector1<potentialtheta3D*> potential_bundle; //list of potentials

ComboPatch(int n) : potential_bundle(vector1<potentialtheta3D*>(n)), safe(true) { p = new int*;}

virtual int num_patches(const int&) = 0; //return the number of patches on particle int
virtual void UpdateIterator(const int&, const int&) = 0; //reassign the pointer to point to something else
virtual void UpdateIteratorSafe(const int&, const int&, int**) =0 ;
virtual int get_total_patches(const int &) = 0;
virtual void which_patch(const int&, const int&, const int &, int&, int & ) = 0; //get patch numbers from particle numbers
virtual void which_particle(const int&, const int&, int&, int&)=0; //get particle numbers from patch numbers
virtual int which_potential(const int&, const int&, const int&, const int& )= 0;
virtual void get_params(const int&,const int&, const int&, double&, double&, double&, double&, double&, double& , double&, double&) = 0;
virtual ComboPatch *clone() const = 0;

virtual void CreateFiles();
//~ComboPatch() {delete p;}


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

    ~SingPatch() {delete i1;}

    int num_patches(const int&) {return 1;}
    void UpdateIterator(const int &i,const int &j) {
    //do nothing (only one patch)
    }
    void UpdateIteratorSafe(const int &i,const int &j, int **q) {
        *q =  i1;
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

    matrix<double> v;

    TetrahedralPatch(double strr, double disss, double angg);

    ~TetrahedralPatch() {
        delete i1;
    }

    int num_patches(const int &i)
    {
        return 4;
    }
    void UpdateIterator(const int &i, const int &j)
    {
        //do nothing (only one patch)
    }
    void UpdateIteratorSafe(const int &i, const int &j, int **q)
    {
        *q = i1;
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

    void CreateFiles();

    TetrahedralPatch *clone() const
    {
        return new TetrahedralPatch(*this);
    }
};


struct TetrahedralWithSingle : ComboPatch {
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

    double angtt;
    double distt;
    double strtt;

    double angts;
    double dists;
    double strts;

    double angss;
    double disss;
    double strss;

    int *i1; // t<-> t
    int *i2;//  t<-> s
    int *i3;//  t<-> s
    
    int nt; //number of tetrahedra
    int ns; //number of singles
    matrix<double> v;

    TetrahedralWithSingle(double, double, double, double, double, double, double, double, double, int nt, int ns);

    ~TetrahedralWithSingle() {
        delete i1;
        delete i2;
        delete i3;
    }

    int num_patches(const int &i)
    {
        if(i < nt ) {
            return 4;
        }
        else if(i >= nt ) { 
            return 1;
        }
        else{
            error("out of bounds");
        }
    }
    void UpdateIterator(const int &i, const int &j)
    {
     if(i < nt && j < nt ) {
         p = &i1;
     }
     else if (i >= nt && j >= nt ) {
         p = &i3;
     }
     else {
         p = &i2;
     }
    }
    void UpdateIteratorSafe(const int &i, const int &j, int **q)
    {
        if (i < nt && j < nt)
        {
            *q = i1;
        }
        else if (i >= nt && j >= nt)
        {
            *q = i3;
        }
        else
        {
            *q = i2;
        }
    }
    int get_total_patches(const int &N) { return 4*nt+ns; }

    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
    {
        if(i < nt && j < nt ) {
        int k1 = potn / 4;
        int k2 = potn % 4;
        // i*4 + k1;
        // j*4 + k2;
        wpi = i * 4 + k1;
        wpj = j * 4 + k2;
        }
        else if(i >= nt && j >= nt) {
            wpi = nt * 4 + (i - nt);
            wpj = nt * 4 + (j - nt);
        }
        else if(i < nt && j >= nt) {
            int k1 = potn % 4;
            wpi = i * 4 + k1;
            wpj = nt * 4 + (j - nt);
        }
        else {
            int k2 =  potn % 4;
            wpi = nt *4 + (i - nt);
            wpj = j*4 + k2;
        }
    }
    void which_particle(const int &wpi, const int &wpj, int &i, int &j)
    {
        if(wpi < 4*nt && wpj < 4*nt) {
            i = wpi / 4;
            j = wpj / 4;
        }
        else if(wpi >= 4*nt && wpj >= 4*nt ) {
            i = wpi - 4*nt + nt;
            j = wpj - 4*nt + nt;
        }
        else if(wpi < 4*nt && wpj >= 4*nt ){
            i = wpi/4;
            j = wpj - 4*nt + nt;
        }
        else{
            i = wpi - 4*nt + nt;
            j = wpj/4;
        }
    }
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
    {
        if (wpi < 4 * nt && wpj < 4 * nt)
        {
            return (wpi % 4) * 4 + (wpj % 4);
        }
        else if (wpi >= 4 * nt && wpj >= 4 * nt)
        {
            return 20;
        }
        else if (wpi < 4 * nt && wpj >= 4 * nt)
        {
            return 16+ wpi % 4;
        }
        else{
            return 16+ wpj % 4;
        }
    }
    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);

    void CreateFiles();

    TetrahedralWithSingle *clone() const
    {
        return new TetrahedralWithSingle(*this);
    }
};

struct TwoTetrahedral : ComboPatch
{
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

    double angtt;
    double distt;
    double strtt;

    double angts;
    double dists;
    double strts;

    double angss;
    double disss;
    double strss;

    int nt;
    int ns;

    int *i1;
    int *i2;
    int *i3;

    TwoTetrahedral(double, double, double, double, double, double, double, double, double, int nt, int ns);

    ~TwoTetrahedral()
    {
        delete i1;
        delete i2;
        delete i3;
    }

    int num_patches(const int &i)
    {
        return 4;
    }
    void UpdateIterator(const int &i, const int &j)
    {
        if (i < nt && j < nt)
        {
            p = &i1;
        }
        else if (i >= nt && j >= nt)
        {
            p = &i3;
        }
        else
        {
            p = &i2;
        }
    }
    void UpdateIteratorSafe(const int &i, const int &j, int **q)
    {
        if (i < nt && j < nt)
        {
            *q = i1;
        }
        else if (i >= nt && j >= nt)
        {
            *q = i3;
        }
        else
        {
            *q = i2;
        }
    }

    int get_total_patches(const int &N) { return 4 * nt + 4*ns; }

    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
    {
        if(i< nt && j < nt) {
            int k1 = potn / 4;
            int k2 = potn % 4;
            // i*4 + k1;
            // j*4 + k2;
            wpi = i * 4 + k1;
            wpj = j * 4 + k2;
        }
        else if(i >= nt && j >= nt) {
            int potn2 = potn - 32;
            int k1 = potn2 / 4;
            int k2 = potn2 % 4;
            // i*4 + k1;
            // j*4 + k2;
            wpi = i * 4 + k1;
            wpj = j * 4 + k2;
        }
        else
        {
            int potn2 = potn - 16;
            int k1 = potn2 / 4;
            int k2 = potn2 % 4;
            // i*4 + k1;
            // j*4 + k2;
            wpi = i * 4 + k1;
            wpj = j * 4 + k2;
        }
    }

    void which_particle(const int &wpi, const int &wpj, int &i, int &j)
    {

            i = wpi / 4;
            j = wpj / 4;
    }
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
    {
        if (wpi < 4 * nt && wpj < 4 * nt)
        {
            return (wpi % 4) * 4 + (wpj % 4);
        }
        else if (wpi >= 4 * nt && wpj >= 4 * nt)
        {
            return 32+(wpi % 4) * 4 + (wpj % 4);
        }
        else
        {
            return 16 + (wpi % 4) * 4 + (wpj % 4);
        }
    }
    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);

    void CreateFiles();

    TwoTetrahedral *clone() const
    {
        return new TwoTetrahedral(*this);
    }
};

struct TwoTetrahedralAndSingle : ComboPatch
{
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

    // double angt1t1;
    // double dist1t1;
    // double strt1t1;

    // double angt1t2;
    // double dist1t2;
    // double strt1t2;

    // double angt1s;
    // double dist1s;
    // double strt1s;

    // double angt2t2;
    // double dist2t2;
    // double strt2t2;

    // double angt2s;
    // double dist2s;
    // double strt2s;

    // double angss;
    // double disss;
    // double strss;

    matrix<double> params;

    int nt;
    int ns;
    int nf;

    int *i1; //p1 p1
    int *i2; //p1 p2
    int *i3; //p1 p3
    int *i4; //p2 p2
    int *i5; //p2 p3
    int *i6; //p3 p3

    matrix<double> v;

    TwoTetrahedralAndSingle(matrix<double>&, int nt, int ns, int nf);

    ~TwoTetrahedralAndSingle()
    {
        delete i1;
        delete i2;
        delete i3;
        delete i4;
        delete i5;
        delete i6;
    }

    int num_patches(const int &i)
    {
        if(i<ns) return 4;
        else return 1;
    }

    inline int mapping_funcion_particles(int i, int j) {
        if (i < nt)
        {
            if (j < nt)
            {
               return 1;
            }
            else if (j < ns)
            {
                return 2;
            }
            else
            {
                return 3;
            }
        }
        else if (i < ns)
        {
            if (j < nt)
            {
                return 2;
            }
            else if (j < ns)
            {
                return 4;
            }
            else
            {
                return 5;
            }
        }
        else
        {
            if (j < nt)
            {
                return 3;
            }
            else if (j < ns)
            {
                return 5;
            }
            else
            {
                return 6;
            }
        }
    }
    void UpdateIterator(const int &i, const int &j)
    {
        int we_are = mapping_funcion_particles(i,j);

        switch(we_are) {
            case 1:
                p = &i1;
                break;

            case 2:
                p = &i2;
                break;

            case 3:
                p = &i3;
                break;

            case 4:
                p = &i4;
                break;

            case 5:
                p = &i5;
                break;

            case 6:
                p = &i6;
                break;
        }
    }
    void UpdateIteratorSafe(const int &i, const int &j, int **q)
    {
        int we_are = mapping_funcion_particles(i, j);

        switch (we_are)
        {
        case 1:
            *q = i1;
            break;

        case 2:
            *q = i2;
            break;

        case 3:
            *q = i3;
            break;

        case 4:
            *q = i4;
            break;

        case 5:
            *q = i5;
            break;

        case 6:
            *q = i6;
            break;
        }
    }

    int get_total_patches(const int &N) { return 4 * nt + 4 * (ns-nt) + nf; }

    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
    {
        int we_are = mapping_funcion_particles(i, j);

        switch (we_are)
        {
        case 1: {
            int k1 = potn / 4;
            int k2 = potn % 4;
            wpi = i * 4 + k1;
            wpj = j * 4 + k2;
            break;
        }
        case 2: {
            int potn2 = potn - 16;
            int k1 = potn2 / 4;
            int k2 = potn2 % 4;
            wpi = i * 4 + k1;
            wpj = j * 4 + k2;
            break;
        }
        case 3: {
            int k1 = (potn - 32) % 4;
            if(i<j) {
            wpi = i * 4 + k1;
            wpj = ns * 4 + (j - ns);
            }
            else{
                wpj = j * 4 + k1;
                wpi = ns * 4 + (i - ns);
            }
            break;
        }
        case 4: {
            int potn2 = potn - 36;
            int k1 = potn2 / 4;
            int k2 = potn2 % 4;
            wpi = i * 4 + k1;
            wpj = j * 4 + k2;
            break;
        }
        case 5: {
            int k1 = (potn - 52) % 4;
            if (i < j)
            {
                wpi = i * 4 + k1;
                wpj = ns * 4 + (j - ns);
            }
            else
            {
                wpj = j * 4 + k1;
                wpi = ns * 4 + (i - ns);
            }
            break;
        }
        case 6: {
            wpi = ns * 4 + (i - ns);
            wpj = ns * 4 + (j - ns);
            break;
        }
        
        }
        // if (i < nt && j < nt)
        // {
        //     int k1 = potn / 4;
        //     int k2 = potn % 4;
        //     // i*4 + k1;
        //     // j*4 + k2;
        //     wpi = i * 4 + k1;
        //     wpj = j * 4 + k2;
        // }
        // else if (i >= nt && j >= nt)
        // {
        //     wpi = nt * 4 + (i - nt);
        //     wpj = nt * 4 + (j - nt);
        // }
        // else if (i < nt && j >= nt)
        // {
        //     int k1 = potn % 4;
        //     wpi = i * 4 + k1;
        //     wpj = nt * 4 + (j - nt);
        // }
        // else
        // {
        //     int k2 = potn % 4;
        //     wpi = nt * 4 + (i - nt);
        //     wpj = j * 4 + k2;
        // }

        // if (i < nt && j < nt)
        // {
        //     int k1 = potn / 4;
        //     int k2 = potn % 4;
        //     // i*4 + k1;
        //     // j*4 + k2;
        //     wpi = i * 4 + k1;
        //     wpj = j * 4 + k2;
        // }
        // else if (i >= nt && j >= nt)
        // {
        //     int potn2 = potn - 32;
        //     int k1 = potn2 / 4;
        //     int k2 = potn2 % 4;
        //     // i*4 + k1;
        //     // j*4 + k2;
        //     wpi = i * 4 + k1;
        //     wpj = j * 4 + k2;
        // }
        // else
        // {
        //     int potn2 = potn - 16;
        //     int k1 = potn2 / 4;
        //     int k2 = potn2 % 4;
        //     // i*4 + k1;
        //     // j*4 + k2;
        //     wpi = i * 4 + k1;
        //     wpj = j * 4 + k2;
        // }
    }

    void which_particle(const int &wpi, const int &wpj, int &i, int &j)
    {
        if (wpi < 4 * ns && wpj < 4 * ns)
        {
            i = wpi / 4;
            j = wpj / 4;
        }
        else if (wpi >= 4 * ns && wpj >= 4 * ns)
        {
            i = wpi - 4 * ns + ns;
            j = wpj - 4 * ns + ns;
        }
        else if (wpi < 4 * ns && wpj >= 4 * ns)
        {
            i = wpi / 4;
            j = wpj - 4 * ns + ns;
        }
        else
        {
            i = wpi - 4 * ns + ns;
            j = wpj / 4;
        }
    }
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
    {
        int we_are = mapping_funcion_particles(i, j);
        switch (we_are)
        {
        case 1:
            return (wpi % 4) * 4 + (wpj % 4);
            break;

        case 2:
            return 16 + (wpi % 4) * 4 + (wpj % 4);
            break;

        case 3:
            if(i<j) return  32 + wpi % 4;
            else return 32 + wpj % 4;
            break;

        case 4:
            return 36 + (wpi % 4) * 4 + (wpj % 4);
            break;

        case 5:
            if (i < j)
                return 52 + wpi % 4;
            else
                return 52 + wpj % 4;
            break;

        case 6:
            return 56;
            break;
        }

        // if (wpi < 4 * nt && wpj < 4 * nt)
        // {
        //     return (wpi % 4) * 4 + (wpj % 4);
        // }
        // else if (wpi >= 4 * nt && wpj >= 4 * nt)
        // {
        //     return 20;
        // }
        // else if (wpi < 4 * nt && wpj >= 4 * nt)
        // {
        //     return 16 + wpi % 4;
        // }
        // else
        // {
        //     return 16 + wpj % 4;
        // }

        // if (wpi < 4 * nt && wpj < 4 * nt)
        // {
        //     return (wpi % 4) * 4 + (wpj % 4);
        // }
        // else if (wpi >= 4 * nt && wpj >= 4 * nt)
        // {
        //     return 0;
        // }
        // else if (wpi < 4 * nt && wpj >= 4 * nt)
        // {
        //     return 16 + wpi % 4;
        // }
        // else
        // {
        //     return 16 + wpj % 4;
        // }

    }
    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);

    void CreateFiles();

    TwoTetrahedralAndSingle *clone() const
    {
        return new TwoTetrahedralAndSingle(*this);
    }
};

#include "combopatch.cpp"
#include "combopatchoutput.cpp"

#endif /* COMBOPATCH_H */
