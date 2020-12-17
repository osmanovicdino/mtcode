#ifndef BINDINGMODELTERNARY_H
#define BINDINGMODELTERNARY_H

#include "BindingModeAbstract.h"
#include "../../Basic/basic.h"


struct BindingModelTernary : AbstractBindingModel
{
    vector1<double> doubr11;
    vector1<double> doubr12;
    vector1<double> doubr22;
    vector1<double> doubr13;
    vector1<double> doubr23;
    vector1<double> doubr33;

    vector1<double> tripr111;
    vector1<double> tripr112;
    vector1<double> tripr113;
    vector1<double> tripr122;
    vector1<double> tripr123;
    vector1<double> tripr133;
    vector1<double> tripr222;
    vector1<double> tripr223;
    vector1<double> tripr233;
    vector1<double> tripr333;

    int div1;
    int div2;

    BindingModelTernary(int, int);

    void setup(double st11, double st22, double st33, double st12, double st13, double st23, double, double, double, double, double, double, double, double, double);

    inline vector1<double> get_drate(int, int);
    inline vector1<double> get_trate(int,int,int);

    void doublet(bool before, int index1, int index2, bool &after);

    void triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &after1, bool &after2, bool &after3);

    void nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters);

    BindingModelTernary *clone() const
    {
        return new BindingModelTernary(*this);
    }
};

#include "BindingModelTernary.cpp"

#endif /* BINDINGMODELTERNARY_H */
