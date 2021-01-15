#ifndef TETRAHEDRALWITHSINGLE_H
#define TETRAHEDRALWITHSINGLE_H

struct TetrahedralWithSingle : ComboPatch
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

    int *i1; // t<-> t
    int *i2; //  t<-> s
    int *i3; //  t<-> s

    int nt; //number of tetrahedra
    int ns; //number of singles
    matrix<double> v;

    TetrahedralWithSingle(double, double, double, double, double, double, double, double, double, int nt, int ns);

    ~TetrahedralWithSingle()
    {
        delete i1;
        delete i2;
        delete i3;
    }

    int num_patches(const int &i);
    void UpdateIterator(const int &i, const int &j);
    void UpdateIteratorSafe(const int &i, const int &j, int **q);
    int get_total_patches(const int &N);

    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj);
    void which_particle(const int &wpi, const int &wpj, int &i, int &j);
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj);

    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);

    void CreateFiles();

    TetrahedralWithSingle *clone() const
    {
        return new TetrahedralWithSingle(*this);
    }
};



#endif /* TETRAHEDRALWITHSINGLE_H */
