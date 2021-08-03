#ifndef BINDINGMODEABSTRACT_CPP
#define BINDINGMODEABSTRACT_CPP

void  AbstractBindingModel::doublet_eq(bool before, int index1, int index2, bool &after, double energ, double &r) {
    //r is a random number which is between 0 and 1
    
    double bef =  ( (int)(before) ) * energ;

//double r =  (double)rand()/(double)(RAND_MAX);

    double aft =  ( (int)(!before) ) * energ;

    double de = aft - bef;

    if(de < 0) {
        after = !before;
    }
    else if(exp(-de) > r) {
        after = !before;
    }
    else{
        after = before;
    }

}

double energy_choice_from(double e1, double e2, double e3 , int cho){
    switch(cho) {
    case 0:
        return e1;
        break;
    case 1:
        return e2;
        break;
    case 2:
        return e3;
        break;
    case 3:
        return 0;
        break;
    default:
        error("error in energy choice");
    }

    error("something very wrong if got here");
    return e1;
}

void set_from(bool &a12, bool &a23, bool &a13, int cho)
{
    switch (cho) {
    case 0:
        a12 = true;
        a23 = false;
        a13 = false;
        break;
    case 1:
        a12 = false;
        a23 = true;
        a13 = false;
        break;
    case 2:
        a12 = false;
        a23 = false;
        a13 = true;
        break;
    case 3:
        a12 = false;
        a23 = false;
        a13 = false;
        break;
    default:
        cout << cho << endl;
        error("error in setting bool");
    }
}

void AbstractBindingModel::triplet_eq(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &after1, bool &after2, bool &after3, double e12, double e23, double e13, double &r)
{

    double energy_before = ((int)b12) * e12 + ((int)b23) * e23 + ((int)b13) * e13;

    //how many possible states depends on the connectivities
    int got;
    double energy_after;
    got = rand() % 4;

    if(c12 == true && c23 == false && c13 == true ) 
    {
        
        while(got == 1) {
            got = rand() % 4;
        }
        energy_after =  energy_choice_from(e12,e23,e13,got);


    }
    else if (c12 == false && c23 == true && c13 == true)
    {
        
        while (got == 0)
        {
            got = rand() % 4;
        }
        energy_after = energy_choice_from(e12, e23, e13, got);
    }
    else if (c12 == true && c23 == true && c13 == false )
    {
        
        while (got == 2)
        {
            got = rand() % 4;
        }
        energy_after = energy_choice_from(e12, e23, e13, got);
    }
    else if (c12 == true && c23 == true && c13 == true)
    {
        energy_after = energy_choice_from(e12, e23, e13, got);
    }
    else {
        error("error in the connectivity");
    }

    double de = energy_after - energy_before;
    //double r = (double)rand()/(double)(RAND_MAX);

    if(de<0) {
        set_from(after1, after2, after3, got);

    }
    else if(exp(-de) > r ) {
        set_from(after1, after2, after3, got);
    }
    else{
        after1= b12;
        after2= b23;
        after3= b13;
    }
}

#endif /* BINDINGMODEABSTRACT_CPP */
