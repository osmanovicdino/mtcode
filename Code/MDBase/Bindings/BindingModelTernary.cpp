#ifndef BINDINGMODELTERNARY_CPP
#define BINDINGMODELTERNARY_CPP

template <typename Q>
BindingModelTernary<Q>::BindingModelTernary(/* int n, int m,  */Q &func2) : doubr11(vector1<double>(4)),
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

{// {
//     div1 = n;
//     div2 = m;
    func = func2;
}

template<typename Q>
void BindingModelTernary<Q>::setup(double st11, double st22, double st33, double st12, double st13, double st23,
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
    double baserates = 0.1;
    double baserates2 = 1.0;
    double baserates3 = 0.001;

    //base rates will be related to the double rates

    double pos_11 = log((1. - st11) / st11);
    double pos_12 = log((1. - st12) / st12);
    double pos_22 = log((1. - st22) / st22);
    double pos_13 = log((1. - st13) / st13);
    double pos_23 = log((1. - st23) / st23);
    double pos_33 = log((1. - st33) / st33);

    

    // set_stable_triple(tripr111, baserates3, baserates3, baserates2, baserates3, baserates2, baserates2, 0.0, 0.0, pos_11, 0.0, pos_11, pos_11);
    // set_stable_triple(tripr112, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_11_12 , assym_11_12 ,pos_11, 0.0 /*assym_12_12*/  , pos_12, pos_12);
    // set_stable_triple(tripr113, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_11_13, assym_11_13, pos_11, 0.0 /*assym_13_13*/, pos_13, pos_13);
    // set_stable_triple(tripr122, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_12_22 , 0.0 /*assym_12_12*/ , pos_12, -assym_12_22  , pos_22,pos_12);
    // set_stable_triple(tripr123, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_12_23  , assym_12_13  , pos_12, assym_23_13  , pos_23, pos_13);
    // set_stable_triple(tripr133, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_13_33  , 0.0 /*assym_13_13*/  , pos_13, -assym_13_33  , pos_33, pos_13);
    // set_stable_triple(tripr222, baserates, baserates, baserates2, baserates, baserates2, baserates2, 0.0, 0.0, pos_22, 0.0, pos_22, pos_22);
    // set_stable_triple(tripr223, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_22_23  , assym_22_23  , pos_22, 0.0 /*assym_23_23*/  , pos_23, pos_23);
    // set_stable_triple(tripr233, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_23_33  , 0.0 /*assym_23_23*/ , pos_23, -assym_23_33  , pos_23, pos_23);
    // set_stable_triple(tripr333, baserates, baserates, baserates2, baserates, baserates2, baserates2, 0.0, 0.0, pos_33, 0.0, pos_33, pos_33);

    double base_sub_11 = baserates * st11;
    double base_sub_12 = baserates * st12;
    double base_sub_13 = baserates * st13;
    double base_sub_22 = baserates * st22;
    double base_sub_23 = baserates * st23;
    double base_sub_33 = baserates * st33;

    set_stable_triple(tripr111, base_sub_11, base_sub_11, baserates2, base_sub_11, baserates2, baserates2, 0.0, 0.0, pos_11, 0.0, pos_11, pos_11);
    set_stable_triple(tripr112, base_sub_12, base_sub_12, baserates2, base_sub_12, baserates2, baserates2, assym_11_12, assym_11_12, pos_11, 0.0 /*assym_12_12*/, pos_12, pos_12);
    set_stable_triple(tripr113, base_sub_13, base_sub_13, baserates2, base_sub_13, baserates2, baserates2, assym_11_13, assym_11_13, pos_11, 0.0 /*assym_13_13*/, pos_13, pos_13);
    set_stable_triple(tripr122, base_sub_22, base_sub_12, baserates2, base_sub_12, baserates2, baserates2, assym_12_22, 0.0 /*assym_12_12*/, pos_12, -assym_12_22, pos_22, pos_12);
    set_stable_triple(tripr123, base_sub_23, base_sub_13, baserates2, base_sub_13, baserates2, baserates2, assym_12_23, assym_12_13, pos_12, assym_23_13, pos_23, pos_13);
    set_stable_triple(tripr133, base_sub_33, base_sub_13, baserates2, base_sub_13, baserates2, baserates2, assym_13_33, 0.0 /*assym_13_13*/, pos_13, -assym_13_33, pos_33, pos_13);
    set_stable_triple(tripr222, base_sub_22, base_sub_22, baserates2, base_sub_22, baserates2, baserates2, 0.0, 0.0, pos_22, 0.0, pos_22, pos_22);
    set_stable_triple(tripr223, base_sub_23, base_sub_23, baserates2, base_sub_23, baserates2, baserates2, assym_22_23, assym_22_23, pos_22, 0.0 /*assym_23_23*/, pos_23, pos_23);
    set_stable_triple(tripr233, base_sub_33, base_sub_23, baserates2, base_sub_23, baserates2, baserates2, assym_23_33, 0.0 /*assym_23_23*/, pos_23, -assym_23_33, pos_23, pos_23);
    set_stable_triple(tripr333, base_sub_33, base_sub_33, baserates2, base_sub_33, baserates2, baserates2, 0.0, 0.0, pos_33, 0.0, pos_33, pos_33);

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

double baserates2(double st11on, double st11off) {
    if(st11on < 1E-12) st11on = 1E-12;
    if(1-st11on < 1E-12 ) st11on =  1-1E-12;
    if(st11off < 1E-12) st11off = 1E-12;
    if(1-st11off < 1E-12 ) st11off =  1-1E-12;

    return 1. / sqrt(SQR(log(1 - st11on)) - 2 * st11off * SQR(log(1 - st11on)) +SQR(st11off) * SQR(log(1 - st11on)));
}

double myexp(double st11on, double st11off) {
    if(st11on < 1.E-12) st11on = 1.E-12;
    if(1.-st11on < 1.E-12 ) st11on =  1.-1.E-12;
    if(st11off < 1E-12) st11off = 1.E-12;
    if(1.-st11off < 1.E-12 ) st11off =  1.-1.E-12;



    return log(1.-st11off);
}

template <typename Q>
void BindingModelTernary<Q>::setup_energy_barrier(
    double st11_on, double st22_on, double st33_on, double st12_on, double st13_on, double st23_on,
    double st11_off, double st22_off, double st33_off, double st12_off, double st13_off, double st23_off,
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

    //in the presence of an energy barrier, it may be that the forward reaction is slow, but also essentially irreversible

    doubr11[0] = 1. - (double)st11_on; //from unbind to unbind
    doubr11[1] = (double)st11_on;      //from unbind to bind
    doubr11[2] = 1. - (double)st11_off; //from to bind to unbind
    doubr11[3] = (double)st11_off;      //from bind to bind

    doubr22[0] = 1.0 - (double)st22_on;
    doubr22[1] = (double)st22_on;
    doubr22[2] = 1.0 - (double)st22_off;
    doubr22[3] = (double)st22_off;

    doubr33[0] = 1.0 - (double)st33_on;
    doubr33[1] = (double)st33_on;
    doubr33[2] = 1.0 - (double)st33_off;
    doubr33[3] = (double)st33_off;

    doubr12[0] = 1.0 - (double)st33_on;
    doubr12[1] = (double)st33_on;
    doubr12[2] = 1.0 - (double)st33_off;
    doubr12[3] = (double)st33_off;

    doubr13[0] = 1.0 - (double)st13_on;
    doubr13[1] = (double)st13_on;
    doubr13[2] = 1.0 - (double)st13_off;
    doubr13[3] = (double)st13_off;

    doubr23[0] = 1.0 - (double)st23_on;
    doubr23[1] = (double)st23_on;
    doubr23[2] = 1.0 - (double)st23_off;
    doubr23[3] = (double)st23_off;

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
    double baserates = 0.01;
    //double baserates2 = 1.0;
    double baserates3 = 0.001;

    //base rates will be related to the double rates

    // double pos_11 = log((1. - st11_on) / st11_off);
    // double pos_12 = log((1. - st12_on) / st12_off);
    // double pos_22 = log((1. - st22_on) / st22_off);
    // double pos_13 = log((1. - st13_on) / st13_off);
    // double pos_23 = log((1. - st23_on) / st23_off);
    // double pos_33 = log((1. - st33_on) / st33_off);

    double pos_11 = myexp(st11_on, st11_off);
    double pos_12 = myexp(st12_on, st12_off);
    double pos_22 = myexp(st22_on, st22_off);
    double pos_13 = myexp(st13_on, st13_off);
    double pos_23 = myexp(st23_on, st23_off);
    double pos_33 = myexp(st33_on, st33_off);

    // set_stable_triple(tripr111, baserates3, baserates3, baserates2, baserates3, baserates2, baserates2, 0.0, 0.0, pos_11, 0.0, pos_11, pos_11);
    // set_stable_triple(tripr112, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_11_12 , assym_11_12 ,pos_11, 0.0 /*assym_12_12*/  , pos_12, pos_12);
    // set_stable_triple(tripr113, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_11_13, assym_11_13, pos_11, 0.0 /*assym_13_13*/, pos_13, pos_13);
    // set_stable_triple(tripr122, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_12_22 , 0.0 /*assym_12_12*/ , pos_12, -assym_12_22  , pos_22,pos_12);
    // set_stable_triple(tripr123, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_12_23  , assym_12_13  , pos_12, assym_23_13  , pos_23, pos_13);
    // set_stable_triple(tripr133, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_13_33  , 0.0 /*assym_13_13*/  , pos_13, -assym_13_33  , pos_33, pos_13);
    // set_stable_triple(tripr222, baserates, baserates, baserates2, baserates, baserates2, baserates2, 0.0, 0.0, pos_22, 0.0, pos_22, pos_22);
    // set_stable_triple(tripr223, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_22_23  , assym_22_23  , pos_22, 0.0 /*assym_23_23*/  , pos_23, pos_23);
    // set_stable_triple(tripr233, baserates, baserates, baserates2, baserates, baserates2, baserates2, assym_23_33  , 0.0 /*assym_23_23*/ , pos_23, -assym_23_33  , pos_23, pos_23);
    // set_stable_triple(tripr333, baserates, baserates, baserates2, baserates, baserates2, baserates2, 0.0, 0.0, pos_33, 0.0, pos_33, pos_33);

    double base_sub_11 = baserates * st11_on;
    double base_sub_12 = baserates * st12_on;
    double base_sub_13 = baserates * st13_on;
    double base_sub_22 = baserates * st22_on;
    double base_sub_23 = baserates * st23_on;
    double base_sub_33 = baserates * st33_on;


    //WAS PREVIOUSLY (1-STOFF)/(1-STON)
    double base_sub_11_off = -baserates2(st11_on, st11_off) * (1 - st11_off) * log(1 - st11_on);
    double base_sub_12_off = -baserates2(st12_on, st12_off) * (1 - st12_off) * log(1 - st12_on);
    double base_sub_13_off = -baserates2(st13_on, st13_off) * (1 - st13_off) * log(1 - st13_on);
    double base_sub_22_off = -baserates2(st22_on, st22_off) * (1 - st22_off) * log(1 - st22_on);
    double base_sub_23_off = -baserates2(st23_on, st23_off) * (1 - st23_off) * log(1 - st23_on);
    double base_sub_33_off = -baserates2(st33_on, st33_off) * (1 - st33_off) * log(1 - st33_on);

  


    set_stable_triple(tripr111, base_sub_11, base_sub_11, base_sub_11_off, base_sub_11, base_sub_11_off, base_sub_11_off, 0.0, 0.0, pos_11, 0.0, pos_11, pos_11);
    set_stable_triple(tripr112, base_sub_12, base_sub_12, base_sub_11_off, base_sub_12, base_sub_12_off, base_sub_12_off, assym_11_12, assym_11_12, pos_11, 0.0 /*assym_12_12*/, pos_12, pos_12);
    set_stable_triple(tripr113, base_sub_13, base_sub_13, base_sub_11_off, base_sub_13, base_sub_13_off, base_sub_13_off, assym_11_13, assym_11_13, pos_11, 0.0 /*assym_13_13*/, pos_13, pos_13);
    set_stable_triple(tripr122, base_sub_22, base_sub_12, base_sub_12_off, base_sub_12, base_sub_22_off, base_sub_12_off, assym_12_22, 0.0 /*assym_12_12*/, pos_12, -assym_12_22, pos_22, pos_12);
    set_stable_triple(tripr123, base_sub_23, base_sub_13, base_sub_12_off, base_sub_13, base_sub_23_off, base_sub_13_off, assym_12_23, assym_12_13, pos_12, assym_23_13, pos_23, pos_13);
    set_stable_triple(tripr133, base_sub_33, base_sub_13, base_sub_13_off, base_sub_13, base_sub_33_off, base_sub_13_off, assym_13_33, 0.0 /*assym_13_13*/, pos_13, -assym_13_33, pos_33, pos_13);
    set_stable_triple(tripr222, base_sub_22, base_sub_22, base_sub_22_off, base_sub_22, base_sub_22_off, base_sub_22_off, 0.0, 0.0, pos_22, 0.0, pos_22, pos_22);
    set_stable_triple(tripr223, base_sub_23, base_sub_23, base_sub_22_off, base_sub_23, base_sub_23_off, base_sub_22_off, assym_22_23, assym_22_23, pos_22, 0.0 /*assym_23_23*/, pos_23, pos_23);
    set_stable_triple(tripr233, base_sub_33, base_sub_23, base_sub_23_off, base_sub_23, base_sub_23_off, base_sub_23_off, assym_23_33, 0.0 /*assym_23_23*/, pos_23, -assym_23_33, pos_23, pos_23);
    set_stable_triple(tripr333, base_sub_33, base_sub_33, base_sub_33_off, base_sub_33, base_sub_33_off, base_sub_33_off, 0.0, 0.0, pos_33, 0.0, pos_33, pos_33);

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

template <typename Q>
inline vector1<double>& BindingModelTernary<Q>::get_drate(int i, int j)
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


template <typename Q>
void BindingModelTernary<Q>::doublet(bool before, int index1, int index2, bool &after)
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
        //int ind1, ind2;

        int ind1 = func(index1);
        int ind2 = func(index2);

        // cout << ind1 << " " << ind2 << endl;

    
        
        /*if (index1 < div1)
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
            ind2 = 3; */

        int indt1,indt2;
        sort_doublet(ind1,ind2,indt1,indt2);

        rtemp = get_drate(indt1,indt2);



        double r1 = rtemp[i1 * 2 + 0]; //rate to unbound
        double r2 = rtemp[i1 * 2 + 1]; //rate to bound

        double totalr = r1 + r2;

        double rando = ((double)rand() / (double)(RAND_MAX));

        double rr = totalr * ((double)rand() / (double)(RAND_MAX));

        if (rr < r1)
        {
            after = false;
        }
        else
        {
            after = true; //becomes bound
        }

        // if(indt1 == 1 && indt2 == 3 && before == false && after == true) {
        //     cout << ind1 << " " << ind2 << endl;
        //     cout << rtemp << endl;
        //     cout << i1 << endl;
        //     pausel();
        // }

}

template <typename Q>
double BindingModelTernary<Q>::calculate_score(int index1, int index2, bool before)
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
    //int ind1, ind2;

    int ind1 = func(index1);
    int ind2 = func(index2);

    // cout << ind1 << " " << ind2 << endl;

    /*if (index1 < div1)
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
            ind2 = 3; */

    int indt1, indt2;
    sort_doublet(ind1, ind2, indt1, indt2);

    rtemp = get_drate(indt1, indt2);

    if(before) {
    return rtemp[i1 * 2 + 0]; //rate to unbound
    }
    else{
    return rtemp[i1 * 2 + 1]; //rate to bound
    }

    // if(indt1 == 1 && indt2 == 3 && before == false && after == true) {
    //     cout << ind1 << " " << ind2 << endl;
    //     cout << rtemp << endl;
    //     cout << i1 << endl;
    //     pausel();
    // }
}

template <typename Q>
inline vector1<double>& BindingModelTernary<Q>::get_trate(int i, int j, int k) { //MUST BE SORTEd
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

template <class Q>
void BindingModelTernary<Q>::print() {

ofstream myfile;
myfile.open("res.csv");

myfile <<= doubr11;
myfile << "\n";

myfile <<= doubr22;
myfile << "\n";

myfile <<= doubr33;
myfile << "\n";

myfile <<= doubr12;
myfile << "\n";

myfile <<= doubr13;
myfile << "\n";

myfile <<= doubr23;
myfile << "\n";



myfile <<= tripr111;
myfile << "\n";

myfile <<= tripr112;
myfile << "\n";

myfile <<= tripr113;
myfile << "\n";

myfile <<= tripr122;
myfile << "\n";

myfile <<= tripr123;
myfile << "\n";

myfile <<= tripr133;
myfile << "\n";

myfile <<= tripr222;
myfile << "\n";

myfile <<= tripr223;
myfile << "\n";

myfile <<= tripr233;
myfile << "\n";

myfile <<= tripr333;
myfile << "\n";

myfile.close();
}

// stringstream global_triplet_analyze8_12;
// stringstream global_triplet_analyze8_23;
// stringstream global_triplet_analyze8_13;
// stringstream global_triplet_analyze9_12;
// stringstream global_triplet_analyze9_23;
// stringstream global_triplet_analyze9_13;
// stringstream global_triplet_analyze10_12;
// stringstream global_triplet_analyze10_23;
// stringstream global_triplet_analyze10_13;
int totalsl = 0;

int fv(int i) {
    if( i < 4000) {
        return i % 4;
    }
    else {
        return 8 + (i-4000)%4;
    }
}

void output_ss_to_file(string file1, stringstream &ss) {
ofstream myfile;
myfile.open(file1.c_str());
myfile << ss.str();
myfile.close();
}

void append_ss_to_file(string file1, stringstream &ss)
{
    ofstream myfile;
    myfile.open(file1.c_str(), std::ios_base::app);
    myfile << ss.str();
    myfile.close();
}

template <typename Q>
void BindingModelTernary<Q>::triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &a12, bool &a23, bool &a13)
{

    vector1<double> rtemp;
    int ind1 = func(index1);
    int ind2 = func(index2);
    int ind3 = func(index3);






    int indt1,indt2,indt3;
    unsigned char o1,o2,o3;
    sort_triplet_and_save_permutation(ind1,ind2,ind3,indt1,indt2,indt3,o1,o2,o3);



    bool tb12,tb23,tb13;
    bool tc12,tc23,tc13;
    save_to_type(b12, b23, b13, o1, o2, o3, tb12, tb23, tb13);
    save_to_type(c12, c23, c13, o1, o2, o3, tc12, tc23, tc13);

    int i1;
    if (tb12)
    { //INDEX 1 and INDEX2 bound
        i1 = 0;
    }
    else if (tb23)
    { //INDEX 2 and INDEX 3 bound
        i1 = 1;
    }
    else if (tb13)
    { //INDEX 1 and INDEX 3 bound
        i1 = 2;
    }
    else
    { //NO BINDINGS
        i1 = 3;
    }

    // string s1 = to_string(indt1);
    // string s2 = to_string(indt2);
    // string s3 = to_string(indt3);

    // // Concatenate both strings
    // string s = s1 + s2 + s3;

    // // Convert the concatenated string
    // // to integer
    // int str = stoi(s);

    rtemp = get_trate(indt1,indt2,indt3);

    double r1 = tc12 * rtemp[i1 * 4 + 0]; //goes to bind 12
    double r2 = tc23 * rtemp(i1 * 4 + 1); //goes to bind 23
    double r3 = tc13 * rtemp[i1 * 4 + 2]; //goes to bind 13
    double r4 = rtemp[i1 * 4 + 3];       //goes to none-bound




    double totalr = r1 + r2 + r3 + r4;

    double rr = totalr * ((double)rand() / (double)(RAND_MAX));
    

    bool ta12,ta23,ta13;


    if (rr < r1)
    {
        ta12 = true;
        ta23 = false;
        ta13 = false;
    }
    else if (rr >= r1 && rr < r1+r2)
    {
        ta12 = false;
        ta23 = true;
        ta13 = false;
    }
    else if (rr >= r1+r2 && rr < r1+r2+r3)
    {
        ta12 = false;
        ta23 = false;
        ta13 = true;
    }
    else
    {
        ta12 = false;
        ta23 = false;
        ta13 = false;
    }

    //unsigned char w1,w2,w3;
    //what_order(o1,o2,o3,w1,w2,w3);

    inverse_save_to_type(ta12, ta23, ta13, o1, o2, o3, a12, a23, a13);

}

template <typename Q>
void BindingModelTernary<Q>::nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters) {
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

            int ind1 = func(index1);
            int ind2 = func(index2);
/*             int ind1, ind2;

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
                ind2 = 3; */

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
