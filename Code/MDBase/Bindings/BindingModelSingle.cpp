BindingModelSingle::BindingModelSingle(double onn, double offf) : on_rate(onn), off_rate(offf) {
    tr = off_rate + on_rate;
}

void BindingModelSingle::doublet(bool before, int index1, int index2, bool &after)
{


    double totalr = tr;

    double rr = totalr * ((double)rand() / (double)(RAND_MAX));

    if (rr < off_rate)
    {
        after = false;
    }
    else
    {
        after = true; //becomes bound
    }


}

void BindingModelSingle::triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &a12, bool &a23, bool &a13)
{

    double r1;
    double r2;
    double r3;
    double r4;
    int i1;
    if (b12)
    { //INDEX 1 and INDEX2 bound
        i1 = 0;
        r1 = c12*on_rate;
        r2 = c23* ( (off_rate / 3.) + 0.001 );
        r3 = c13*( (off_rate / 3.) + 0.001 );
        r4 = (off_rate / 3.) + 0.001;
    }
    else if (b23)
    { //INDEX 2 and INDEX 3 bound
        i1 = 1;
        r1 = c12 *  ( (off_rate / 3.) + 0.001 );
        r2 = c23 * (on_rate);
        r3 = c13 * ( (off_rate / 3.) + 0.001 );
        r4 = ( (off_rate / 3.) + 0.001 );
    }
    else if (b13)
    { //INDEX 1 and INDEX 3 bound
        i1 = 2;
        r1 = c12 *( (off_rate / 3.) + 0.001 );
        r2 = c23 *( (off_rate / 3.) + 0.001 );
        r3 = c13 * on_rate;
        r4 = (off_rate / 3.) + 0.001;
    }
    else
    { //NO BINDINGS
        i1 = 3;
        r1 = c12 * on_rate/3.;
        r2 = c23 * on_rate/3.;
        r3 = c13 * on_rate/3.;
        r4 = off_rate;
    }


    // double r1 = c12 * triprates(ing, i1 * 4 + 0); //goes to bind 12
    // double r2 = c23 * triprates(ing, i1 * 4 + 1); //goes to bind 23
    // double r3 = c13 * triprates(ing, i1 * 4 + 2); //goes to bind 13
    // double r4 = triprates(ing, i1 * 4 + 3);       //goes to none-bound

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

void BindingModelSingle::nlet(const vector1<bool> &befores, const vector<mdpair> &indices, const vector<vector1<bool>> &possibles, vector1<bool> &afters)
{
    int bb = 0;

    for(int i = 0  ; i < befores.getsize() ; i++) {
        bb += (int)befores.gpcons(i);
    }

    vector1<double> possible_rates(possibles.size());

    for(int i = 0 ; i < possibles.size() ; i++) {
        
        int ab = 0;
        for (int j = 0; j < befores.getsize(); j++)
        {
            ab += (int)possibles[i].gpcons(j);
        }

        if(ab > bb || bb == 0)
        {
            possible_rates[i] = SQR(ab - bb) * on_rate; //transition to more boundedness
         } 
        else if(bb < ab){
             possible_rates[i] = (bb - ab)*off_rate; //less boundedness
        }
        else {
            if(befores == possibles[i] ) {
                possible_rates[i] =  on_rate; //if the same state
            }
            else{
                possible_rates[i] = off_rate; //if different state
            }
        }
    }

    double totalrate = 0.0;
    vector1<double> well_defined_rates(possible_rates.getsize());
    for(int i = 0  ; i < possible_rates.getsize() ; i++) {
        totalrate += possible_rates[i];
        double vertexrate = 0.0;
        for(int k = 0 ; k < i+1 ; k++) {
        vertexrate += possible_rates[k];
        }
        well_defined_rates[i] = vertexrate;
    }

    double rr = totalrate*(double)rand()/(double)(RAND_MAX);

    // for(int i = 0  ; i < well_defined_rates.getsize() ; i++) {
    // cout << well_defined_rates[i] << " ";
    // cout << possibles[i] << endl;
    // }

    // cout << totalrate << endl;
    // cout << rr << endl;
    


    int whichto;
    if(rr < well_defined_rates[0]) {
        whichto = 0;
    }
    else{
        for(int i = 1 ; i < possible_rates.getsize() ; i++) {
            if(rr > well_defined_rates[i-1] && rr < well_defined_rates[i] ) {
                whichto = i;
                break;
            }
        }
    }
    


    afters = possibles[whichto];
    //error("nlets are too hard without simplifications");
}
