#ifndef BINDINGMODELSINGLE_H
#define BINDINGMODELSINGLE_H

#include "BindingModeAbstract.h"

struct BindingModelSingle : AbstractBindingModel { //CAUTION this is too expensive for large systems

//for every configuration, there is a probability to go to something different

double on_rate;
double off_rate;

double tr;


BindingModelSingle(double, double);

void doublet(bool before, int index1, int index2, bool &after);

void triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &after1, bool &after2, bool &after3);

void nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters);

BindingModelSingle *clone() const {
    return new BindingModelSingle(*this);
};

};


#include "BindingModelSingle.cpp"

#endif /* BINDINGMODELSINGLE */
