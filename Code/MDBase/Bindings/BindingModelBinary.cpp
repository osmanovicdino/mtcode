#ifndef BINDINGMODELBINARY_CPP
#define BINDINGMODELBINARY_CPP

BindingModelBinary::BindingModelBinary(int divv) : doubr11(vector1<double>(4)), doubr12(vector1<double>(4)), doubr22(vector1<double>(4)), tripr111(vector1<double>(4 * 4)) ,tripr112(vector1<double>(4 * 4)), tripr122(vector1<double>(4)), tripr222(vector1<double>(4))
{
    div = divv;
}

BindingModelBinary::setup_equilibrium() {

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
   /*  int bb = 0;

    for (int i = 0; i < befores.getsize(); i++)
    {
        bb += (int)befores.gpcons(i);
    }

    vector1<double> possible_rates(possibles.size());

    for (int i = 0; i < possibles.size(); i++)
    {

        int ab = 0;
        for (int j = 0; j < befores.getsize(); j++)
        {
            ab += (int)possibles[i].gpcons(j);
        }


        if (ab > bb)
            possible_rates[i] = (ab - bb) * on_rate; //transition to more boundedness
        else if (bb < ab)
            possible_rates[i] = (bb - ab) * off_rate; //less boundedness
        else
        {
            if (befores == possibles[i])
            {
                possible_rates[i] = on_rate; //if the same state
            }
            else
            {
                possible_rates[i] = off_rate; //if different state
            }
        }
    } */
}

#endif /* BINARYMODELBINARY_CPP */
