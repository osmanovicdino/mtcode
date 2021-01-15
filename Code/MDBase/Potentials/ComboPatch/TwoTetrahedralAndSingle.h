#ifndef TWOTETRAHEDRALANDSINGLE_H
#define TWOTETRAHEDRALANDSINGLE_H

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

    TwoTetrahedralAndSingle(matrix<double> &, int nt, int ns, int nf);

    ~TwoTetrahedralAndSingle()
    {
        delete i1;
        delete i2;
        delete i3;
        delete i4;
        delete i5;
        delete i6;
    }

    int num_patches(const int &i);
    inline int mapping_funcion_particles(int i, int j);
    void UpdateIterator(const int &i, const int &j);
    void UpdateIteratorSafe(const int &i, const int &j, int **q);

    int get_total_patches(const int &N);

    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj);

    void which_particle(const int &wpi, const int &wpj, int &i, int &j);
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj);
    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);

    void CreateFiles();

    TwoTetrahedralAndSingle *clone() const
    {
        return new TwoTetrahedralAndSingle(*this);
    }
};

#include "TwoTetrahedralAndSingle.cpp"

#endif /* TWOTETRAHEDRALANDSINGLE_H */
