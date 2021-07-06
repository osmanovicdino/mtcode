#ifndef BINDINGMODELABSTRACT_H
#define BINDINGMODELABSTRACT_H

#include "../../DataStructures/mdpair.h"
#include "../../DataStructures/vector1.h"

struct AbstractBindingModel {


//index1 and index2 are within critical distance and can bind or not
virtual void doublet(bool before, int index1, int index2, bool &after) = 0;

virtual double calculate_score(int index1, int index2, bool b) = 0;

virtual void triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &after1, bool &after2, bool &after3) = 0;

virtual void nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters) = 0;



virtual AbstractBindingModel *clone() const = 0;

};

void set_stable_triple(vector1<double> &t111,
                       double a, double b, double c, double d, double e, double f,
                       double as, double bs, double cs, double ds, double es, double fs)
{
    // we need a way to quantify the asymmetry
    //


    t111[0] = 1;             //fromi,jtpi,j
    t111[1] = a * exp(as);   //fromi,jtoj,k
    t111[2] = b * exp(bs);   //fromi,jtoi,k
    t111[3] = c * exp(cs);   //fromi,jtonothing
    t111[4] = a * exp(-as);  //fromj,ktoi,j
    t111[5] = 1;             //fromj,ktoj,k
    t111[6] = d * exp(ds);   //fromj,ktoi,k
    t111[7] = e * exp(es);   //fromj,ktonothing
    t111[8] = b * exp(-bs);  //fromi,ktoi,j
    t111[9] = d * exp(-ds);  //fromi,ktoj,k
    t111[10] = 1;            //fromi,ktoi,k
    t111[11] = f * exp(fs);  //fromi,ktonothing
    t111[12] = c * exp(-cs); //fromnothingtoi,j
    t111[13] = e * exp(-es); //fromnothingtoj,k
    t111[14] = f * exp(-fs); //fromnothingtoi,k
    t111[15] = 1.;           //fromnothingtonothing
}

#endif