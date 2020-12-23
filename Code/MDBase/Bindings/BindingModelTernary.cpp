#ifndef BINDINGMODELTERNARY_CPP
#define BINDINGMODELTERNARY_CPP

BindingModelTernary::BindingModelTernary(int n, int m) : doubr11(vector1<double>(4)),
                                                         doubr12(vector1<double>(4)),
                                                         doubr22(vector1<double>(4)),
                                                         doubr13(vector1<double>(4)),
                                                         doubr23(vector1<double>(4)),
                                                         doubr33(vector1<double>(4)),
                                                         tripr111(vector1<double>(16)),
                                                         tripr112(vector1<double>(16)),
                                                         tripr113(vector1<double>(16)),
                                                         tripr122(vector1<double>(16)),
                                                         tripr123(vector1<double>(16)),
                                                         tripr133(vector1<double>(16)),
                                                         tripr222(vector1<double>(16)),
                                                         tripr223(vector1<double>(16)),
                                                         tripr233(vector1<double>(16)),
                                                         tripr333(vector1<double>(16))
{
    div1 = n;
    div2 = m;
}

void BindingModelTernary::setup(double st11, double st22, double st33, double st12, double st13, double st23,
                                double assym_11_12,
                                double assym_11_13,
                                double assym_12_22,
                                double assym_12_23,
                                double assym_12_13,
                                double assym_23_13,
                                double assym_13_33,
                                double assym_22_23,
                                double assym_23_33

)
{

    doubr11[0] = 1.-(double)st11; //from unbind to unbind
    doubr11[1] = (double)st11;   //from unbind to bind
    doubr11[2] = 1.-(double)st11; //from to bind to unbind
    doubr11[3] = (double)st11; //from bind to bind


    doubr22[0] = 1.0-(double)st22;
    doubr22[1] = (double)st22;
    doubr22[2] = 1.0-(double)st22;
    doubr22[3] = (double)st22;

    doubr33[0] = 1.0-(double)st33;
    doubr33[1] = (double)st33;
    doubr33[2] = 1.0-(double)st33;
    doubr33[3] = (double)st33;  

    doubr12[0] = 1.0-(double)st33;
    doubr12[1] = (double)st33;
    doubr12[2] = 1.0-(double)st33;
    doubr12[3] = (double)st33;  

    doubr13[0] = 1.0-(double)st13;
    doubr13[1] = (double)st13;
    doubr13[2] = 1.0-(double)st13;
    doubr13[3] = (double)st13;  

    doubr23[0] = 1.0-(double)st23;
    doubr23[1] = (double)st23;
    doubr23[2] = 1.0-(double)st23;
    doubr23[3] = (double)st23;

    double baserates = 0.1;

    //if st11 is 1, it means it binds strongly, if it is zero it means it unbinds strongly
    //we can represent this assymetry with 1-2*st11, so that  the rate of going strongly to unbound states matches the doublet
    //when all particles are the same, there are no assymetries in the bindings between them, only between different particles
    //there we replace all the factors corresponding to transitions between two molecules bound of the same type with 0

    //for a general form
    //set_stable_triple(triprabc, baserates, baserates, baserates, baserates, baserates, baserates,a,b,c,d,e,f);

    //the parameters c,e,f refer to the asymmetries of (a,b)->nothing, (b,c)->nothing and (a,c)->nothing

    //if this value is positive it means that there is a high rate for bound states to go to unbound state
    //if this value is negative it means vice versa.

    //We are then left with three other parameters, a,b,d

    // a is the asymmetry between (a,b)<->(b,c)
    // b is the asymmetry between (a,b)<->(a,c)
    // d is the asymmetry between (b,c)<->(a,c)

    // double sub_assym_12 = 1.0; //assymetry of (1,1)->(1,2)
    // double sub_assym_21 = 1.0; //assymetry of (2,2) -> (1,2)

    // double sub_assym_13 = 1.0;
    // double assym_11_12;
    // double assym_11_13;
    // double assym_12_22;

    // double assym_12_23;
    // double assym_12_13;
    // double assym_23_13;

    // double assym_13_33;

    // double assym_22_23;
    // double assym_23_33;

    //only keep allowed transitions

    double baserates2 = 1.0;
    double baserates3 = 0.001;
    //base rates will be related to the double rates

    double pos_11 = log((1. - st11) / st11);
    double pos_12 = log((1. - st12) / st12);
    double pos_22 = log((1. - st22) / st22);
    double pos_13 = log((1. - st13) / st13);
    double pos_23 = log((1. - st23) / st23);
    double pos_33 = log((1. - st33) / st33);

    set_stable_triple(tripr111, baserates3, baserates3, baserates2, baserates3, baserates2, baserates2, 0.0, 0.0, pos_11, 0.0, pos_11, pos_11);
    set_stable_triple(tripr112, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_11_12 , assym_11_12 ,pos_11, 0.0 /*assym_12_12*/  , pos_12, pos_12);
    set_stable_triple(tripr113, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_11_13, assym_11_13, pos_11, 0.0 /*assym_13_13*/, pos_13, pos_13);
    set_stable_triple(tripr122, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_12_22 , 0.0 /*assym_12_12*/ , pos_12, -assym_12_22  , pos_22,pos_12);
    set_stable_triple(tripr123, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_12_23  , assym_12_13  , pos_12, assym_23_13  , pos_23, pos_13);
    set_stable_triple(tripr133, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_13_33  , 0.0 /*assym_13_13*/  , pos_13, -assym_13_33  , pos_33, pos_13);
    set_stable_triple(tripr222, baserates, baserates, baserates2, baserates, baserates2, baserates2, 0.0, 0.0, pos_22, 0.0, pos_22, pos_22);
    set_stable_triple(tripr223, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_22_23  , assym_22_23  , pos_22, 0.0 /*assym_23_23*/  , pos_23, pos_23);
    set_stable_triple(tripr233, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_23_33  , 0.0 /*assym_23_23*/ , pos_23, -assym_23_33  , pos_23, pos_23);
    set_stable_triple(tripr333, baserates, baserates, baserates2, baserates, baserates2, baserates2, 0.0, 0.0, pos_33, 0.0, pos_33, pos_33);

    /*

    t111[0] = 0.998;  //fromi,jtpi,j
    t111[1] = 0.001;  //fromi,jtoj,k
    t111[2] = 0.001;  //fromi,jtoi,k
    t111[3] = 0.000;  //fromi,jtonothing
    t111[4] = 0.001;  //fromj,ktoi,j
    t111[5] = 0.998;  //fromj,ktoj,k
    t111[6] = 0.001;  //fromj,ktoi,k
    t111[7] = 0.000;  //fromj,ktonothing
    t111[8] = 0.001;  //fromi,ktoi,j
    t111[9] = 0.001;  //fromi,ktoj,k
    t111[10] = 0.998; //fromi,ktoi,k
    t111[11] = 0.000; //fromi,ktonothing
    t111[12] = 0.33;  //fromnothingtoi,j
    t111[13] = 0.33;  //fromnothingtoj,k
    t111[14] = 0.33;  //fromnothingtoi,k
    t111[15] = 0.01;  //fromnothingtonothing

    //particle 2 does not interact

    t112[0] = 0.001;  //fromi,jtpi,j
    t112[1] = 0.499;  //fromi,jtoj,k
    t112[2] = 0.499;  //fromi,jtoi,k
    t112[3] = 0.000;  //fromi,jtonothing
    t112[4] = 0.000;  //fromj,ktoi,j
    t112[5] = 0.998;  //fromj,ktoj,k
    t112[6] = 0.100;  //fromj,ktoi,k
    t112[7] = 0.001;  //fromj,ktonothing
    t112[8] = 0.000;  //fromi,ktoi,j
    t112[9] = 0.100;  //fromi,ktoj,k
    t112[10] = 0.998; //fromi,ktoi,k
    t112[11] = 0.001; //fromi,ktonothing
    t112[12] = 0.001; //fromnothingtoi,j
    t112[13] = 0.499; //fromnothingtoj,k
    t112[14] = 0.499; //fromnothingtoi,k
    t112[15] = 0.000; //fromnothingtonothing

    t122[0] = 0.998;  //fromi,jtpi,j
    t122[1] = 0.000;  //fromi,jtoj,k
    t122[2] = 0.100;  //fromi,jtoi,k
    t122[3] = 0.001;  //fromi,jtonothing
    t122[4] = 0.499;  //fromj,ktoi,j
    t122[5] = 0.000;  //fromj,ktoj,k
    t122[6] = 0.499;  //fromj,ktoi,k
    t122[7] = 0.001;  //fromj,ktonothing
    t122[8] = 0.100;  //fromi,ktoi,j
    t122[9] = 0.000;  //fromi,ktoj,k
    t122[10] = 0.998; //fromi,ktoi,k
    t122[11] = 0.001; //fromi,ktonothing
    t122[12] = 0.499; //fromnothingtoi,j
    t122[13] = 0.000; //fromnothingtoj,k
    t122[14] = 0.499; //fromnothingtoi,k
    t122[15] = 0.001; //fromnothingtonothing

    t222[0] = 0.000;  //fromi,jtpi,j
    t222[1] = 0.000;  //fromi,jtoj,k
    t222[2] = 0.000;  //fromi,jtoi,k
    t222[3] = 1.0;    //fromi,jtonothing
    t222[4] = 0.000;  //fromj,ktoi,j
    t222[5] = 0.000;  //fromj,ktoj,k
    t222[6] = 0.000;  //fromj,ktoi,k
    t222[7] = 1.000;  //fromj,ktonothing
    t222[8] = 0.000;  //fromi,ktoi,j
    t222[9] = 0.000;  //fromi,ktoj,k
    t222[10] = 0.000; //fromi,ktoi,k
    t222[11] = 1.000; //fromi,ktonothing
    t222[12] = 0.000; //fromnothingtoi,j
    t222[13] = 0.000; //fromnothingtoj,k
    t222[14] = 0.000; //fromnothingtoi,k
    t222[15] = 1.000; //fromnothingtonothing

    */
}

inline vector1<double> BindingModelTernary::get_drate(int i, int j)
{
    if(i==j) {
        if(i==1) {
            return doubr11;
        }
        else if(i==2) {
            return doubr22;
        }
        else{
            return doubr33;
        }
    }
    else{
        if(j==2) {
        return doubr12;
        }
        else{
            if(i==2) {
                return doubr23;
            }
            else{
                return doubr13;
            }
        }
    }
    
}

void BindingModelTernary::doublet(bool before, int index1, int index2, bool &after)
{


        int i1;
        if (before == false)
        {
            i1 = 0;
        }
        else
        {
            i1 = 1;
        }

        vector1<double> rtemp;
        int ind1, ind2;

        if (index1 < div1)
            ind1 = 1;
        else if (index1 < div2)
            ind1 = 2;
        else
            ind1 = 3;

        if (index2 < div1)
            ind2 = 1;
        else if (index2 < div2)
            ind2 = 2;
        else
            ind2 = 3;

        int indt1,indt2;
        sort_doublet(ind1,ind2,indt1,indt2);

        rtemp = get_drate(indt1,indt2);

        double r1 = rtemp[i1 * 2 + 0]; //rate to unbound
        double r2 = rtemp[i1 * 2 + 1]; //rate to bound

        double totalr = r1 + r2;

        double rr = totalr * ((double)rand() / (double)(RAND_MAX));

        if (rr < r1)
        {
            after = false;
        }
        else
        {
            after = true; //becomes bound
        }


}


inline vector1<double> BindingModelTernary::get_trate(int i, int j, int k) { //MUST BE SORTEd
    if(i==1) {
        if(j==1) {
            if(k==1) {
                return tripr111;
            }
            else if(k==2) {
                return tripr112;
            }
            else{
                return tripr113;
            }
        }
        else if(j==2){
            if(k==2) {
                return tripr122;
            }
            else{
                return tripr123;
            }
        }
        else{
            return tripr133;
        }
    }
    else if(i==2) {
        if(j==2) {
            if(k==2) {
                return tripr222;
            }
            else{
                return tripr223;
            }
        }
        else{
            return tripr233;
        }

    }
    else{
        return tripr333;
    }
}

void BindingModelTernary::triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &a12, bool &a23, bool &a13)
{
    int i1;
    if (b12)
    { //INDEX 1 and INDEX2 bound
        i1 = 0;
    }
    else if (b23)
    { //INDEX 2 and INDEX 3 bound
        i1 = 1;
    }
    else if (b13)
    { //INDEX 1 and INDEX 3 bound
        i1 = 2;
    }
    else
    { //NO BINDINGS
        i1 = 3;
    }
    vector1<double> rtemp;
    int ind1,ind2,ind3;

    if(index1 < div1) 
        ind1 = 1;
    else if(index1 < div2)
        ind1 = 2;
    else 
        ind1 = 3;

    if (index2 < div1)
        ind2 = 1;
    else if (index2 < div2)
        ind2 = 2;
    else
        ind2 = 3;

    if (index3 < div1)
        ind3 = 1;
    else if (index3 < div2)
        ind3 = 2;
    else
        ind3 = 3;

    int indt1,indt2,indt3;
    sort_triplet(ind1,ind2,ind3,indt1,indt2,indt3);

    rtemp = get_trate(indt1,indt2,indt3);

    

    double r1 = c12 * rtemp[i1 * 4 + 0]; //goes to bind 12
    double r2 = c23 * rtemp(i1 * 4 + 1); //goes to bind 23
    double r3 = c13 * rtemp[i1 * 4 + 2]; //goes to bind 13
    double r4 = rtemp[i1 * 4 + 3];       //goes to none-bound

    double totalr = r1 + r2 + r3 + r4;

    double rr = totalr * ((double)rand() / (double)(RAND_MAX));

    if (rr < r1)
    {
        a12 = true;
        a23 = false;
        a13 = false;
    }
    else if (rr >= r1 && rr < r2)
    {
        a12 = false;
        a23 = true;
        a13 = false;
    }
    else if (rr >= r2 && rr < r3)
    {
        a12 = false;
        a23 = false;
        a13 = true;
    }
    else
    {
        a12 = false;
        a23 = false;
        a13 = false;
    }
}

void BindingModelTernary::nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters) {
    // {
    //     cout << befores << endl;
    //     pausel();
    
    int bb = 0;
    int nb = befores.getsize();


    for (int i = 0; i < nb; i++)
    {
        bb += (int)befores.gpcons(i);
    }

    int nst = possibles.size(); //number of states to

    vector1<double> possible_rates(nst);

    for (int i = 0; i < nst; i++)
    { //for all the number of possible states

        for(int j = 0 ; j < nb ; j++) {
            vector1<double> rtemp;
            int index1 = indices[j].a;
            int index2 = indices[j].b;

            int ind1, ind2;

            if (index1 < div1)
                ind1 = 1;
            else if (index1 < div2)
                ind1 = 2;
            else
                ind1 = 3;

            if (index2 < div1)
                ind2 = 1;
            else if (index2 < div2)
                ind2 = 2;
            else
                ind2 = 3;

            int indt1, indt2;
            sort_doublet(ind1, ind2, indt1, indt2);

            rtemp = get_drate(indt1, indt2);

            int i1;
            if (befores.gpcons(j) == false)
            {
                i1 = 0;
            }
            else
            {
                i1 = 1;
            }

            int i2 = (int)possibles[i].gpcons(j);

            possible_rates[i] += rtemp[i1*2+i2];
            //possible_rates[j] += rtemp();

        }

    }

    double totalrate = 0.0;
    vector1<double> well_defined_rates(nst);
    for (int i = 0; i < nst; i++)
    {
        totalrate += possible_rates[i];
        double vertexrate = 0.0;
        for (int k = 0; k < i + 1; k++)
        {
            vertexrate += possible_rates[k];
        }
        well_defined_rates[i] = vertexrate;
    }

    double rr = totalrate * (double)rand() / (double)(RAND_MAX);



    int whichto;
    if (rr < well_defined_rates[0])
    {
        whichto = 0;
    }
    else
    {
        for (int i = 1; i < nst; i++)
        {
            if (rr > well_defined_rates[i - 1] && rr < well_defined_rates[i])
            {
                whichto = i;
                break;
            }
        }
    }

    afters = possibles[whichto];
    

   //afters = befores;
}

#endif /* BINDINGMODELTERNARY_CPP */
