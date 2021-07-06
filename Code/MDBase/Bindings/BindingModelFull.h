#ifndef BINDINGMODELFULL_H
#define BINDINGMODELFULL_H

#include "BindingModeAbstract.h"

struct BindingModelFull : AbstractBindingModel { //CAUTION this is too expensive for large systems

//for every configuration, there is a probability to go to something different
int N;
matrix<double> doubrates;
matrix<double> triprates;

BindingModelFull(int);

void doublet(bool before, int index1, int index2, bool &after);

double calculate_score(int index1, int index2, bool are_b);

void triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &after1, bool &after2, bool &after3);

void nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters);

BindingModelFull *clone() const {
    return new BindingModelFull(*this);
};

};


#include "BindingModelFull.cpp"

#endif /* BINDINGMODELFULL_H */
