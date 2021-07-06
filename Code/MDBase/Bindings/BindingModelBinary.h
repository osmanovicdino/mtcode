#ifndef BINDINGMODELBINARY_H
#define BINDINGMODELBINARY_H

#include "BindingModeAbstract.h"

struct BindingModelBinary : AbstractBindingModel
{ //CAUTION this is too expensive for large systems

    //for every configuration, there is a probability to go to something different
    vector1<double> doubr11;
    vector1<double> doubr12;
    vector1<double> doubr22;

    vector1<double> tripr111;
    vector1<double> tripr112;
    vector1<double> tripr122;
    vector1<double> tripr222;

    int div;

    BindingModelBinary(int);

    void setup_equilibrium();

    void doublet(bool before, int index1, int index2, bool &after);

    double calculate_score(int index1, int index2, bool b);

    void triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &after1, bool &after2, bool &after3);

    void nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters);

    BindingModelBinary *clone() const
    {
        return new BindingModelBinary(*this);
    }
};

#include "BindingModelBinary.cpp"

#endif /* BINDINGMODELBINARY_H */
