#ifndef TETRAHEDRALWITHBIVALENT_H
#define TETRAHEDRALWITHBIVALENT_H

struct TetrahedralWithBivalent : ComboPatch
{
    // double nx1 = sqrt(8. / 9.);
    // double ny1 = 0.;
    // double nz1 = -1. / 3.;

    // double nx2 = -sqrt(2. / 9.);
    // double ny2 = sqrt(2. / 3.);
    // double nz2 = -1. / 3.;

    // double nx3 = -sqrt(2. / 9.);
    // double ny3 = -sqrt(2. / 3.);
    // double nz3 = -1. / 3.;

    // double nx4 = 0;
    // double ny4 = 0;
    // double nz4 = 1.;

    // double nx5 = -1.0;
    // double ny5 = 0.0;
    // double nz5 = 0.0;

    // double nx6 = 1.0;
    // double ny6 = 0.0;
    // double nz6 = 0.0;

    int *i1; // t<-> t
    int *i2; //  t<-> s
    int *i3; //  t<-> s

    int nt; //number of tetrahedra
    int nb; //number of singles
    matrix<double> v;
    matrix<double> v2;

    matrix<double> params2;

    matrix<double> defaultv() const {
        matrix<double> v(4,3);
        v(0, 0) = sqrt(8. / 9.);
        v(0, 1) = 0.;
        v(0, 2) = -1. / 3.;

        v(1, 0) = -sqrt(2. / 9.);
        v(1, 1) = sqrt(2. / 3.);
        v(1, 2) = -1. / 3.;

        v(2, 0) = -sqrt(2. / 9.);
        v(2, 1) = -sqrt(2. / 3.);
        v(2, 2) = -1. / 3.;

        v(3, 0) = 0;
        v(3, 1) = 0;
        v(3, 2) = 1.;

        return v;
    }

    matrix<double> defaultv2() const
    {
        matrix<double> v2(2, 3);
        v2(0, 0) = -1.;
        v2(0, 1) = 0.;
        v2(0, 2) = 0.;

        v2(1, 0) = 1.;
        v2(1, 1) = 0.;
        v2(1, 2) = 0.;

        return v2;
    }

    TetrahedralWithBivalent();
    TetrahedralWithBivalent(matrix<double> &,int,int);
    TetrahedralWithBivalent(matrix<double> &,int,int, const matrix<double>&,const matrix<double>&);

    ~TetrahedralWithBivalent()
    {
        delete i1;
        delete i2;
        delete i3;
    }

    void change_nt(const int &ntt) { nt = ntt;}

    int num_patches(const int &i);
    void UpdateIterator(const int &i, const int &j);
    void UpdateIteratorSafe(const int &i, const int &j, int **q);
    int get_total_patches(const int &N);

    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj);
    void which_particle(const int &wpi, const int &wpj, int &i, int &j);
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj);

    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);

    void CreateFiles();

    TetrahedralWithBivalent *clone() const
    {
        return new TetrahedralWithBivalent(*this);
    }

};

#include "TetrahedralWithBivalent.cpp"

#endif /* TETRAHEDRALWITHBIVALENT_H */
