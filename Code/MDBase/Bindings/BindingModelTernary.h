#ifndef BINDINGMODELTERNARY_H
#define BINDINGMODELTERNARY_H

#include "BindingModeAbstract.h"
#include "../../Basic/basic.h"


struct SortingFunctionUniform {
int div1;
int div2;

int operator()(const int &index1) {
    if (index1 < div1)
        return 1;
    else if (index1 < div2)
        return 2;
    else
        return 3;
}


};

struct SortingFunctionNonUniform {
    int div1;
    int div2;
    int np;
    int operator()(const int &index1) {
        if (index1 < div1)
            return 1;
        else if (index1 < div2)
            return 2;
        else {
           // cout << index1-div2 << endl;
            if((index1-div2)%np == np-1) return 3;
            else return 1;
        }
    }
};


template <typename Q> //Modified to include the fact that we shall be able to bin into types a,b,c according to an arbitrary rule passed by template
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

    // int div1;
    // int div2;

    Q func; /* this is the function that takes the input argument and outputs what chemical species it is
    we do not do thiss with runtime polymorphism as this function will need to be called so many times
    as to make the difference between the two non-negligible 

    The only feature this function needs is operator()(int) that returns another int 1 2 or 3
    */

    BindingModelTernary(Q&);

    void setup(double st11, double st22, double st33, double st12, double st13, double st23, double, double, double, double, double, double, double, double, double);

    void setup_energy_barrier(double st11, double st22, double st33, double st12, double st13, double st23, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

    vector1<double>& get_drate(int, int);
    vector1<double>& get_trate(int, int, int);

    void doublet(bool before, int index1, int index2, bool &after);

    double calculate_score(int index1, int index2, bool are_b);

    void triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &after1, bool &after2, bool &after3);

    void nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters);


    void print();

    BindingModelTernary *clone() const
    {
        return new BindingModelTernary(*this);
    }
};

#include "BindingModelTernary.cpp"

#endif /* BINDINGMODELTERNARY_H */
