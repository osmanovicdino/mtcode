#ifndef BINDINGMODELFULL_H
#define BINDINGMODELFULL_H

#include "BindingModeAbstract.h"

struct BindingModelFull : AbstractBindingModel { //CAUTION this is too expensive for large systems

//for every configuration, there is a probability to go to something different
int N;
matrix<double> doubrates;
matrix<double> triprates;

void doublet(bool before, int index1, int index2, bool &after);

void triplet(bool fc, bool bef1, bool bef2, bool bef3, int index1, int index2, int index3, bool &after1, bool after2, bool &after3);

void nlet(vector1<bool> &befores, vector1<int> indices, vector1<bool> &afters) ;

};

#endif /* BINDINGMODELFULL_H */
