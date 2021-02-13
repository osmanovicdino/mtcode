#ifndef BINDINGMODELBINARY_CPP
#define BINDINGMODELBINARY_CPP

BindingModelBinary::BindingModelBinary(int divv) : doubr11(vector1<double>(4)), doubr12(vector1<double>(4)), doubr22(vector1<double>(4)), tripr111(vector1<double>(4 * 4)) ,tripr112(vector1<double>(4 * 4)), tripr122(vector1<double>(4)), tripr222(vector1<double>(4))
{
    div = divv;
}

void BindingModelBinary::setup_equilibrium() {

    //particles type 1 and 2 can bind to each other;

    vector1<double> r11(4);
    r11[0] = 0.0;
    r11[1] = 1.0;
    r11[2] = 0.0;
    r11[3] = 1.0;

    vector1<double> r12(r11);
    vector1<double> r22(r11);
    // r22[0] = 1.0;
    // r22[1] = 0.0;
    // r22[2] = 1.0;
    // r22[3] = 0.0;

    vector1<double> t111(16);
    vector1<double> t112(16);
    vector1<double> t122(16);
    vector1<double> t222(16);

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
    t112[0] = 0.998;  //fromi,jtpi,j
    t112[1] = 0.001;  //fromi,jtoj,k
    t112[2] = 0.001;  //fromi,jtoi,k
    t112[3] = 0.000;  //fromi,jtonothing
    t112[4] = 0.001;  //fromj,ktoi,j
    t112[5] = 0.998;  //fromj,ktoj,k
    t112[6] = 0.001;  //fromj,ktoi,k
    t112[7] = 0.000;  //fromj,ktonothing
    t112[8] = 0.001;  //fromi,ktoi,j
    t112[9] = 0.001;  //fromi,ktoj,k
    t112[10] = 0.998; //fromi,ktoi,k
    t112[11] = 0.000; //fromi,ktonothing
    t112[12] = 0.33;  //fromnothingtoi,j
    t112[13] = 0.33;  //fromnothingtoj,k
    t112[14] = 0.33;  //fromnothingtoi,k
    t112[15] = 0.01;  //fromnothingtonothing

    t122[0] = 0.998;  //fromi,jtpi,j
    t122[1] = 0.001;  //fromi,jtoj,k
    t122[2] = 0.001;  //fromi,jtoi,k
    t122[3] = 0.000;  //fromi,jtonothing
    t122[4] = 0.001;  //fromj,ktoi,j
    t122[5] = 0.998;  //fromj,ktoj,k
    t122[6] = 0.001;  //fromj,ktoi,k
    t122[7] = 0.000;  //fromj,ktonothing
    t122[8] = 0.001;  //fromi,ktoi,j
    t122[9] = 0.001;  //fromi,ktoj,k
    t122[10] = 0.998; //fromi,ktoi,k
    t122[11] = 0.000; //fromi,ktonothing
    t122[12] = 0.33;  //fromnothingtoi,j
    t122[13] = 0.33;  //fromnothingtoj,k
    t122[14] = 0.33;  //fromnothingtoi,k
    t122[15] = 0.01;  //fromnothingtonothing

    t222[0] = 0.998;  //fromi,jtpi,j
    t222[1] = 0.001;  //fromi,jtoj,k
    t222[2] = 0.001;  //fromi,jtoi,k
    t222[3] = 0.000;  //fromi,jtonothing
    t222[4] = 0.001;  //fromj,ktoi,j
    t222[5] = 0.998;  //fromj,ktoj,k
    t222[6] = 0.001;  //fromj,ktoi,k
    t222[7] = 0.000;  //fromj,ktonothing
    t222[8] = 0.001;  //fromi,ktoi,j
    t222[9] = 0.001;  //fromi,ktoj,k
    t222[10] = 0.998; //fromi,ktoi,k
    t222[11] = 0.000; //fromi,ktonothing
    t222[12] = 0.33;  //fromnothingtoi,j
    t222[13] = 0.33;  //fromnothingtoj,k
    t222[14] = 0.33;  //fromnothingtoi,k
    t222[15] = 0.01;  //fromnothingtonothing
   
    doubr11 = r11;
    doubr12 = r12;
    doubr22 = r22;

    tripr111 = t111;
    tripr112 = t112;
    tripr122 = t122;
    tripr222 = t222;
}

void BindingModelBinary::doublet(bool before, int index1, int index2, bool &after)
{

    if(index1 < div && index2 < div) {
        int i1;
        if (before == false)
        {
            i1 = 0;
        }
        else
        {
            i1 = 1;
        }
        double r1 =doubr11[ i1 * 2 + 0]; //rate to unbound
        double r2 =doubr11[ i1 * 2 + 1]; //rate to bound

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
    else if(index1 < div && index2 >= div) {
        int i1;
        if (before == false)
        {
            i1 = 0;
        }
        else
        {
            i1 = 1;
        }
        double r1 = doubr12[i1 * 2 + 0]; //rate to unbound
        double r2 = doubr12[i1 * 2 + 1]; //rate to bound

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
    else{
        int i1;
        if (before == false)
        {
            i1 = 0;
        }
        else
        {
            i1 = 1;
        }
        double r1 = doubr22[i1 * 2 + 0]; //rate to unbound
        double r2 = doubr22[i1 * 2 + 1]; //rate to bound

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

}

void BindingModelBinary::triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &a12, bool &a23, bool &a13) {

    if(index1 < div && index2 < div && index3 < div) {
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

        double r1 = c12 * tripr111[i1 * 4 + 0]; //goes to bind 12
        double r2 = c23 * tripr111(i1 * 4 + 1);  //goes to bind 23
        double r3 = c13 * tripr111[i1 * 4 + 2];  //goes to bind 13
        double r4 = tripr111[i1 * 4 + 3];        //goes to none-bound



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
    else if(index1 < div && index2 < div && index3 >= div) {
        //cout << "triple replace event called" << endl;
        
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


        // cout << index1 << endl;
        // cout << index2 << endl;
        // cout << index3 << endl;

        // cout << b12 << " " << b23 << " " << b13 << endl;

        // cout << c12 << " " << c23 << " " << c13 << endl;
        

        double r1 = c12 * tripr112[i1 * 4 + 0]; //goes to bind 12
        double r2 = c23 * tripr112(i1 * 4 + 1); //goes to bind 23
        double r3 = c13 * tripr112[i1 * 4 + 2]; //goes to bind 13
        double r4 = tripr112[i1 * 4 + 3];       //goes to none-bound

        double totalr = r1 + r2 + r3 + r4;

        // cout << totalr << endl;
        // cout << r1 << " " << r2 << " " << r3 << " " << r4 << endl;


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

        // cout << a12 << " " << a23 << " " << a13 << endl;
        // cout << "triple event" << endl;
        // pausel();
    }
    else if(index1 < div && index2 >= div && index3 >= div) {
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

        double r1 = c12 * tripr122[i1 * 4 + 0]; //goes to bind 12
        double r2 = c23 * tripr122(i1 * 4 + 1); //goes to bind 23
        double r3 = c13 * tripr122[i1 * 4 + 2]; //goes to bind 13
        double r4 = tripr111[i1 * 4 + 3];       //goes to none-bound


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
    else {
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

        double r1 = c12 * tripr222[i1 * 4 + 0]; //goes to bind 12
        double r2 = c23 * tripr222(i1 * 4 + 1); //goes to bind 23
        double r3 = c13 * tripr222[i1 * 4 + 2]; //goes to bind 13
        double r4 = tripr111[i1 * 4 + 3];       //goes to none-bound

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


}

void BindingModelBinary::nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool> > &possibles, vector1<bool> &afters)
{
    // int bb = 0;
    // int nb = befores.getsize();

    // for (int i = 0; i < nb; i++)
    // {
    //     bb += (int)befores.gpcons(i);
    // }

    // int nst = possibles.size(); //number of states to

    // vector1<double> possible_rates(nst);

    // for (int i = 0; i < nst; i++)
    // { //for all the number of possible states

    //     for (int j = 0; j < nb; j++)
    //     {
    //         vector1<double> rtemp;
    //         int index1 = indices[j].a;
    //         int index2 = indices[j].b;

    //         int ind1, ind2;

    //         if (index1 < div)
    //             ind1 = 1;
    //         else
    //             ind1 = 2;

    //         if (index2 < div)
    //             ind2 = 1;
    //         else
    //             ind2 = 2;
            
    //         int indt1, indt2;
    //         sort_doublet(ind1, ind2, indt1, indt2);

    //         //rtemp = get_drate(indt1, indt2);

    //         int i1;
    //         if (befores.gpcons(j) == false)
    //         {
    //             i1 = 0;
    //         }
    //         else
    //         {
    //             i1 = 1;
    //         }

    //         int i2 = (int)possibles[i].gpcons(j);
    //         if(indt1 == 1 && indt2 == 1) {
    //             possible_rates[i] += doubr11[i1 * 2 + i2];
    //         }
    //         else if (indt1 == 1 && indt2 == 2) {
    //             possible_rates[i] += doubr12[i1 * 2 + i2];
    //         }
    //         else {
    //             possible_rates[i] += doubr22[i1 * 2 + i2];
    //         }
    //         //possible_rates[j] += rtemp();
    //     }
    // }

    // double totalrate = 0.0;
    // vector1<double> well_defined_rates(nst);
    // for (int i = 0; i < nst; i++)
    // {
    //     totalrate += possible_rates[i];
    //     double vertexrate = 0.0;
    //     for (int k = 0; k < i + 1; k++)
    //     {
    //         vertexrate += possible_rates[k];
    //     }
    //     well_defined_rates[i] = vertexrate;
    // }

    // double rr = totalrate * (double)rand() / (double)(RAND_MAX);

    // int whichto;
    // if (rr < well_defined_rates[0])
    // {
    //     whichto = 0;
    // }
    // else
    // {
    //     for (int i = 1; i < nst; i++)
    //     {
    //         if (rr > well_defined_rates[i - 1] && rr < well_defined_rates[i])
    //         {
    //             whichto = i;
    //             break;
    //         }
    //     }
    // }

    // afters = possibles[whichto];
}

#endif /* BINARYMODELBINARY_CPP */
