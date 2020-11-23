#ifndef BINDINGMODELABSTRACT_H
#define BINDINGMODELABSTRACT_H

#include "../../DataStructures/mdpair.h"

struct AbstractBindingModel {


//index1 and index2 are within critical distance and can bind or not
virtual void doublet(bool before, int index1, int index2, bool &after) = 0;

virtual void triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &after1, bool &after2, bool &after3) = 0;

virtual void nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters) = 0;

virtual AbstractBindingModel *clone() const = 0;

};

#endif