#ifndef BINDINGMODELFULL_CPP
#define BINDINGMODELFULL_CPP

BindingModelFull::BindingModelFull(int n) : N(n), doubrates(matrix<double>(n * n,2)), triprates(matrix<double>(n * n * n,16))
{

}

void BindingModelFull::doublet(bool before, int index1, int index2, bool &after) {
    //the doublet model takes as a function the original state and indices, and produces a final state
    //INDEX1 and INDEX2 are SORTED

    int i1;
    if(before == false) {
        i1 = 0;
    }
    else {
        i1 = 1;
    }

    int ing = index1*N+index2;
    // int rate_index = index1*N+index2;
    // if(before = true) {
    // double rate = doubrates_off[rate_index];
    
    // double r1 = (double)rand()/(double)(RAND_MAX)
    // if(rate > r1 ) {
    //     after = true;
    // }
    // else {
    //     after = false;
    // }
    // }
    // else {
    // double rate =  doubrates_on[rate_index];

    // double r1 = (double)rand()(double)(RAND_MAX);

    // if(rate > r1 ) {
    //     after = true;
    // }
    // else{
    //     after = false;
    // }

    // }
}


//FULL CHARACERITZATION WOULD NEED US TO DEFINE 16 matrices for every possible combo here
void BindingModelFull::triplet(bool b12, bool b23, bool b13, bool c12, bool c23, bool c13, int index1, int index2, int index3, bool &a12, bool &a23, bool &a13)
{

    //MAKE SURE THAT THE INDICES ARE SORTED IN ASCENDING ORDER
    //take the original state
    //i1,i2,i3 are the original indices
    //b13,b13,b23 are the initial bindings. Note, at most only one of these can be true
    //c13,c23,c23 are the connections, defining whether the graph is fully connected or not

    //from these functions we get the rates to a new (consistent) state

    //DETERMINE THE APPROPRIATE RATE MATRIX:
    
    //INITIAL STATE;
    int i1;
    if(b12) { //INDEX 1 and INDEX2 bound
        i1=0;
    }
    else if(b23) { //INDEX 2 and INDEX 3 bound
        i1 = 1;
    }
    else if(b13) { //INDEX 1 and INDEX 3 bound
        i1 = 2;
    }
    else{ //NO BINDINGS
        i1 = 3;
    }

    //this characterizes the initial state

    int ing =  N*N*index1 + N*index2 + index3;

    //WE DEFINE THIS IN THE SAME WAY WE WOULD A MATRIX

    //rates for this given state

    double r1 = c12*triprates(ing, i1 * 4 + 0); //goes to bind 12
    double r2 = c23*triprates(ing, i1 * 4 + 1); //goes to bind 23
    double r3 = c13*triprates(ing, i1 * 4 + 2); //goes to bind 13
    double r4 = triprates(ing, i1 * 4 + 3); //goes to none-bound

    double totalr = r1+r2+r3+r4;

    double rr = totalr*(double)rand()/(double)(RAND_MAX);

    if(rr < r1) {
        a12 = true;
        a23 = false;
        a13 = false;
    }
    else if(rr>=r1 && rr <r2) {
        a12 = false;
        a23 = true;
        a13 = false;
    }
    else if(rr>=r2 && rr < r3 ) {
        a12 = false;
        a23 = false;
        a13 = true;
    }
    else{
        a12 = false;
        a23 = false;
        a13 = false;
    }


    //this gives us the rates for each initial state:
    

    // else 
    // {
    //     //NOT FULLY CONNECTED
    //     //depening on which one is missing
    //     if(!c12) { //if c12 is the connection that is lacking
    //         a12 = false;

    //     }
    //     else if(!c13) {

    //     }
    //     else if(!c23) {

    //     }
    //     else{

    //     }
    // }

}

void BindingModelFull::nlet(vector1<bool> &befores, vector1<int> indices, vector1<bool> &afters)
{
    error("nlets are too hard without simplifications");

}

#endif /* BINDINGMODELFULL_CPP */
