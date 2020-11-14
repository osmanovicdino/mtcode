Connected components example

<pre><code>
matrix<int> adj(200,4);
vector1<int> lens(200);

int tlens[200] ={0,2,1,2,1,0,1,0,1,0,0,0,1,0,0,1,1,2,1,1,1,0,1,2,0,1,0,0,1,2,0,2,0,2,1,0,0,0,2,4,2,0,0,1,0,2,1,1,0,1,0,0,0,1,1,0,0,1,1,1,1,0,1,1,0,1,2,0,1,1,3,0,0,0,3,0,0,3,1,0,1,1,0,1,2,1,2,0,1,1,0,1,0,1,1,0,1,1,2,1,0,2,0,0,0,0,2,0,1,0,0,0,1,2,0,2,0,0,1,2,2,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,2,0,1,0,0,0,0,1,2,1,0,0,1,1,0,0,0,1,0,0,0,1,0,1,0,2,1,1,0,2,0,0,0,0,1,2,0,1,0,0,1,1,1,0,1,1,2,1,1,1,0,1,1,0,1,0,1,2,1,2,2,1,0,4};

 int tadj[200][4] = {{0,0,0,0},{96,132,0,0},{199,0,0,0},{153,162,0,0},{165,0,0,0},{0,0,0,0},{170,0,0,0},{0,0,0,0},{134,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{45,0,0,0},{0,0,0,0},{0,0,0,0},{74,0,0,0},{159,0,0,0},{77,188,0,0},{120,0,0,0},{77,0,0,0},{101,0,0,0},{0,0,0,0},{122,0,0,0},{136,185,0,0},{0,0,0,0},{74,0,0,0},{0,0,0,0},{0,0,0,0},{40,0,0,0},{173,187,0,0},{0,0,0,0},{39,65,0,0},{0,0,0,0},{176,190,0,0},{106,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{60,199,0,0},{31,49,91,99},{28,89,0,0},{0,0,0,0},{0,0,0,0},{69,0,0,0},{0,0,0,0},{12,196,0,0},{163,0,0,0},{59,0,0,0},{0,0,0,0},{39,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{70,0,0,0},{98,0,0,0},{0,0,0,0},{0,0,0,0},{98,0,0,0},{113,0,0,0},{47,0,0,0},{38,0,0,0},{0,0,0,0},{84,0,0,0},{171,0,0,0},{0,0,0,0},{31,0,0,0},{181,199,0,0},{0,0,0,0},{138,0,0,0},{43,0,0,0},{53,119,182,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{15,25,106,0},{0,0,0,0},{0,0,0,0},{17,19,84,0},{144,0,0,0},{0,0,0,0},{101,0,0,0},{149,0,0,0},{0,0,0,0},{86,0,0,0},{62,77,0,0},{194,0,0,0},{83,192,0,0},{0,0,0,0},{144,0,0,0},{40,0,0,0},{0,0,0,0},{39,0,0,0},{0,0,0,0},{193,0,0,0},{197,0,0,0},{0,0,0,0},{1,0,0,0},{145,0,0,0},{54,57,0,0},{39,0,0,0},{0,0,0,0},{20,80,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{34,74,0,0},{0,0,0,0},{171,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{157,0,0,0},{58,195,0,0},{0,0,0,0},{178,184,0,0},{0,0,0,0},{0,0,0,0},{182,0,0,0},{70,193,0,0},{18,165,0,0},{0,0,0,0},{22,0,0,0},{177,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{1,0,0,0},{0,0,0,0},{8,0,0,0},{0,0,0,0},{23,143,0,0},{0,0,0,0},{68,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{136,0,0,0},{78,88,0,0},{97,0,0,0},{0,0,0,0},{0,0,0,0},{195,0,0,0},{81,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{112,0,0,0},{0,0,0,0},{16,0,0,0},{0,0,0,0},{183,196,0,0},{3,0,0,0},{46,0,0,0},{0,0,0,0},{4,120,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{6,0,0,0},{63,108,0,0},{0,0,0,0},{29,0,0,0},{0,0,0,0},{0,0,0,0},{33,0,0,0},{123,0,0,0},{115,0,0,0},{0,0,0,0},{199,0,0,0},{66,0,0,0},{70,118,0,0},{161,0,0,0},{115,0,0,0},{23,0,0,0},{0,0,0,0},{29,0,0,0},{17,0,0,0},{0,0,0,0},{33,0,0,0},{0,0,0,0},{86,0,0,0},{93,119,0,0},{85,0,0,0},{113,148,0,0},{45,161,0,0},{94,0,0,0},{0,0,0,0},{2,38,66,180}};

for (int i = 0; i < 200; i++)
{
    lens[i] = tlens[i];
}

for (int i = 0; i < 200; i++)
for(int j = 0 ; j < 4 ; j++) {
    {
    adj(i,j) = tadj[i][j];
    }
}
//FUNC CC

vector1<int> indexes(200);

vector1<int> nbins = ConnectedComponents(adj,lens,indexes);

for(int i = 0 ; i < nbins.getsize()-1 ; i++) {
    for(int j = nbins[i] ; j < nbins[i+1] ; j++) {
        cout << indexes[j] << " ";
    }
    cout << endl;

}
</code></pre>

This will print out all the connected components for the defined adjacency graph

Apparently working model with single binding condition (equilibrium)

<pre><code>


srand (time(NULL));

int n = 4;

Condensate A(5.0, n); //Create a condensate


BindingModelFull b(n);


matrix<double> mat1(n*n,4);

for(int i = 0  ; i < n ; i++) {
    for(int j = i+1  ; j < n ; j++) {
        int index1 = i*n+j;
        mat1(index1,0) = 0.0; //unbound->unbound
        mat1(index1,1) = 1.0; //unbound->bound
        mat1(index1,2) = 0.0; //bound -> unbound
        mat1(index1,3) = 1.0; //bound -> bound
    }
}





b.doubrates = mat1;


matrix<double> tm(4,4);

tm(0,0) = 0.9;
tm(0,1) = 0.05;
tm(0,2) = 0.05;
tm(0,3) = 0.0;
tm(1,0) = 0.05;
tm(1,1) = 0.9;
tm(1,2) = 0.05;
tm(1,3) = 0.0;
tm(2,0) = 0.05;
tm(2,1) = 0.05;
tm(2,2) = 0.9;
tm(2,3) = 0.0;
tm(3,0) = 0.0;
tm(3,1) = 0.0;
tm(3,2) = 0.0;
tm(3,3) = 0.0;

matrix<double> triprates(n*n*n,16);
for(int i  = 0  ; i < n ; i++) {
    for(int j = i+1 ; j < n ; j++) {
        for(int k = j+1 ; k < n ; k++) {
          //  cout << i << " " << j << " " << k << endl;
            int indx = i*n*n+j*n+k;
            //cout << indx << endl;

            triprates(indx, 0)  = tm(0, 0); // = 0.9;
            triprates(indx, 1) = tm(0, 1);  // = 0.05;
            triprates(indx, 2) = tm(0, 2);  // = 0.05;
            triprates(indx, 3) = tm(0, 3); // = 0.0;
            triprates(indx, 4) = tm(1, 0); // = 0.05;
            triprates(indx, 5) = tm(1, 1); // = 0.9;
            triprates(indx, 6) = tm(1, 2); // = 0.05;
            triprates(indx, 7) = tm(1, 3); // = 0.0;
            triprates(indx, 8) = tm(2, 0); // = 0.05;
            triprates(indx, 9) = tm(2, 1); // = 0.05;
            triprates(indx, 10) = tm(2, 2); // = 0.9;
            triprates(indx, 11) = tm(2, 3); // = 0.0;
            triprates(indx, 12) = tm(3, 0); // = 0.0;
            triprates(indx, 13) = tm(3, 1); // = 0.0;
            triprates(indx, 14) = tm(3, 2); // = 0.0;
            triprates(indx, 15) = tm(3, 3); //= 0.0;
        }
    }
}


b.triprates =  triprates;

SingPatch c(100.0,2.,pi/3.);

A.setBindingModel(b);

A.setpots(c);

A.obj->setkT(1.0);


A.run_singlebond( 1000000,  1000);

</code></pre>

for a binding circulation:

<pre><code>

int n = 4;

Condensate A(5.0, n); //Create a condensate


BindingModelFull b(n);


matrix<double> matb(n*n,4);

for(int i = 0  ; i < n ; i++) {
    for(int j = i+1  ; j < n ; j++) {
        int index1 = i*n+j;
        if(i == 0) { //we define particle 0 as the promiscuous particle (can bind to anything)
            matb(index1,0) = 0.0; //unbound->unbound
            matb(index1,1) = 1.0; //unbound->bound
            matb(index1,2) = 0.0; //bound -> unbound
            matb(index1,3) = 1.0; //bound -> bound
        }
        else{ //no other particles can bind
            matb(index1, 0) = 1.0; //unbound->unbound
            matb(index1, 1) = 0.0; //unbound->bound
            matb(index1, 2) = 1.0; //bound -> unbound
            matb(index1, 3) = 0.0; //bound -> bound
        }
    }
}





b.doubrates = matb;



matrix<double> mat1(4, 4);
matrix<double> mat2(4, 4);
matrix<double> mat3(4, 4);
matrix<double> mat4(4, 4);


//i ==0  j ==1 k == 2
mat1(0,0) = 0.001;
mat1(0,1) = 0.001;
mat1(0,2) = 0.998;
mat1(0,3) = 0.0;
mat1(1,0) = 0.001;
mat1(1,1) = 0.001;
mat1(1,2) = 0.098;
mat1(1,3) = 0.00;
mat1(2,0) = 0.001;
mat1(2,1) = 0.001;
mat1(2,2) = 0.998;
mat1(2,3) = 0.0;
mat1(3,0) = 0.001;
mat1(3,1) = 0.001;
mat1(3,2) = 0.998;
mat1(3,3) = 0.0;

//i ==0  j ==1 k == 3
mat2(0, 0) = 0.998;
mat2(0, 1) = 0.001;
mat2(0, 2) = 0.001;
mat2(0, 3) = 0.0;
mat2(1, 0) = 0.998;
mat2(1, 1) = 0.001;
mat2(1, 2) = 0.001;
mat2(1, 3) = 0.00;
mat2(2, 0) = 0.998;
mat2(2, 1) = 0.001;
mat2(2, 2) = 0.001;
mat2(2, 3) = 0.0;
mat2(3, 0) = 0.998;
mat2(3, 1) = 0.001;
mat2(3, 2) = 0.001;
mat2(3, 3) = 0.0;


//i ==0  j ==2 k == 3
mat3(0, 0) = 0.001;
mat3(0, 1) = 0.001;
mat3(0, 2) = 0.998;
mat3(0, 3) = 0.0;
mat3(1, 0) = 0.001;
mat3(1, 1) = 0.001;
mat3(1, 2) = 0.998;
mat3(1, 3) = 0.00;
mat3(2, 0) = 0.001;
mat3(2, 1) = 0.001;
mat3(2, 2) = 0.998;
mat3(2, 3) = 0.0;
mat3(3, 0) = 0.001;
mat3(3, 1) = 0.001;
mat3(3, 2) = 0.998;
mat3(3, 3) = 0.0;

//i ==1  j ==2 k == 3
mat4(0, 0) = 0.0;
mat4(0, 1) = 0.0;
mat4(0, 2) = 0.0;
mat4(0, 3) = 1.0;
mat4(1, 0) = 0.0;
mat4(1, 1) = 0.0;
mat4(1, 2) = 0.0;
mat4(1, 3) = 1.0;
mat4(2, 0) = 0.0;
mat4(2, 1) = 0.0;
mat4(2, 2) = 0.0;
mat4(2, 3) = 1.0;
mat4(3, 0) = 0.0;
mat4(3, 1) = 0.0;
mat4(3, 2) = 0.0;
mat4(3, 3) = 1.0;

matrix<double> triprates(n*n*n,16);
for(int i  = 0  ; i < n ; i++) {
    for(int j = i+1 ; j < n ; j++) {
        for(int k = j+1 ; k < n ; k++) {
          //  cout << i << " " << j << " " << k << endl;
            int indx = i*n*n+j*n+k;
            //cout << indx << endl;
            if(i == 0 && j ==1 && k == 2) {
                triprates(indx, 0) =  mat1(0, 0); //from i,j tp i,j
                triprates(indx, 1) =  mat1(0, 1);  //from i,j to j,k
                triprates(indx, 2) =  mat1(0, 2);  //from i,j to i,k
                triprates(indx, 3) =  mat1(0, 3);  //from i,j to nothing
                triprates(indx, 4) =  mat1(1, 0);  //from j,k to i,j
                triprates(indx, 5) =  mat1(1, 1);  //from j,k to j,k
                triprates(indx, 6) =  mat1(1, 2);  //from j,k to i,k
                triprates(indx, 7) =  mat1(1, 3);  //from j,k to nothing
                triprates(indx, 8) =  mat1(2, 0);  //from i,k to i,j
                triprates(indx, 9) =  mat1(2, 1);  //from i,k to j,k
                triprates(indx, 10) = mat1(2, 2); //from i,k to i,k
                triprates(indx, 11) = mat1(2, 3); //from i,k to nothing
                triprates(indx, 12) = mat1(3, 0); //from nothing to i,j
                triprates(indx, 13) = mat1(3, 1); //from nothing to j,k
                triprates(indx, 14) = mat1(3, 2); //from nothing to i,k
                triprates(indx, 15) = mat1(3, 3); //from nothing to nothing
            }
            else if(i==0 && j==1 && k ==3) 
            {
                triprates(indx, 0) = mat2(0, 0);  //from i,j tp i,j
                triprates(indx, 1) = mat2(0, 1);  //from i,j to j,k
                triprates(indx, 2) = mat2(0, 2);  //from i,j to i,k
                triprates(indx, 3) = mat2(0, 3);  //from i,j to nothing
                triprates(indx, 4) = mat2(1, 0);  //from j,k to i,j
                triprates(indx, 5) = mat2(1, 1);  //from j,k to j,k
                triprates(indx, 6) = mat2(1, 2);  //from j,k to i,k
                triprates(indx, 7) = mat2(1, 3);  //from j,k to nothing
                triprates(indx, 8) = mat2(2, 0);  //from i,k to i,j
                triprates(indx, 9) = mat2(2, 1);  //from i,k to j,k
                triprates(indx, 10) = mat2(2, 2); //from i,k to i,k
                triprates(indx, 11) = mat2(2, 3); //from i,k to nothing
                triprates(indx, 12) = mat2(3, 0); //from nothing to i,j
                triprates(indx, 13) = mat2(3, 1); //from nothing to j,k
                triprates(indx, 14) = mat2(3, 2); //from nothing to i,k
                triprates(indx, 15) = mat2(3, 3); //from nothing to nothing
            }
            else if (i == 0 && j == 2 && k == 3)
            {
                triprates(indx, 0) = mat3(0, 0);  //from i,j tp i,j
                triprates(indx, 1) = mat3(0, 1);  //from i,j to j,k
                triprates(indx, 2) = mat3(0, 2);  //from i,j to i,k
                triprates(indx, 3) = mat3(0, 3);  //from i,j to nothing
                triprates(indx, 4) = mat3(1, 0);  //from j,k to i,j
                triprates(indx, 5) = mat3(1, 1);  //from j,k to j,k
                triprates(indx, 6) = mat3(1, 2);  //from j,k to i,k
                triprates(indx, 7) = mat3(1, 3);  //from j,k to nothing
                triprates(indx, 8) = mat3(2, 0);  //from i,k to i,j
                triprates(indx, 9) = mat3(2, 1);  //from i,k to j,k
                triprates(indx, 10) = mat3(2, 2); //from i,k to i,k
                triprates(indx, 11) = mat3(2, 3); //from i,k to nothing
                triprates(indx, 12) = mat3(3, 0); //from nothing to i,j
                triprates(indx, 13) = mat3(3, 1); //from nothing to j,k
                triprates(indx, 14) = mat3(3, 2); //from nothing to i,k
                triprates(indx, 15) = mat3(3, 3); //from nothing to nothing
            }
            else if (i == 1 && j == 2 && k == 3)
            {
                triprates(indx, 0) = mat4(0, 0);  //from i,j tp i,j
                triprates(indx, 1) = mat4(0, 1);  //from i,j to j,k
                triprates(indx, 2) = mat4(0, 2);  //from i,j to i,k
                triprates(indx, 3) = mat4(0, 3);  //from i,j to nothing
                triprates(indx, 4) = mat4(1, 0);  //from j,k to i,j
                triprates(indx, 5) = mat4(1, 1);  //from j,k to j,k
                triprates(indx, 6) = mat4(1, 2);  //from j,k to i,k
                triprates(indx, 7) = mat4(1, 3);  //from j,k to nothing
                triprates(indx, 8) = mat4(2, 0);  //from i,k to i,j
                triprates(indx, 9) = mat4(2, 1);  //from i,k to j,k
                triprates(indx, 10) = mat4(2, 2); //from i,k to i,k
                triprates(indx, 11) = mat4(2, 3); //from i,k to nothing
                triprates(indx, 12) = mat4(3, 0); //from nothing to i,j
                triprates(indx, 13) = mat4(3, 1); //from nothing to j,k
                triprates(indx, 14) = mat4(3, 2); //from nothing to i,k
                triprates(indx, 15) = mat4(3, 3); //from nothing to nothing
            }
            else{
                //
            }
        }
    }
}



b.triprates =  triprates;

SingPatch c(100.0,2.,pi/3.);

A.setBindingModel(b);

A.setpots(c);

A.obj->setkT(1.0);


A.run_singlebond( 100000,  100);

</code></pre>