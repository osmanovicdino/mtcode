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

void BindingModelSingle::nlet(vector1<bool> &befores, vector1<int> indices, vector1<bool> &afters)
{
    error("nlets are too hard without simplifications");
}
