https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/



# Connected components example

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

# Apparently working model with single binding condition (equilibrium)

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

# Working Binding Model:

<pre><code>
int n = 125;

Condensate A(21.5443, n);


//for one hundeed twenty five 
double basex =  10.0;
double basey =  10.0;
double basez =  10.0;

double dx =1.0;

int iter  = 0;
matrix<double> initialpos(n,3);

for(int i = 0  ; i < 5 ; i++) {
    for(int j =  0 ; j < 5 ; j ++) {
        for(int k = 0 ; k < 5 ; k++ ) {
            initialpos(iter, 0) = basex + i * dx;
            initialpos(iter, 1) = basey + j * dx;
            initialpos(iter, 2) = basez + k * dx;
            iter++;
        }
    }
}

A.obj->setdat(initialpos);

TetrahedralPatch c(30.0, 1.4, pi / 4.);

BindingModelSingle b(1.0,0.0);

A.setBindingModel(b);

A.setpots(c);


for (double kT = 1.0; kT > 0.49; kT -= 0.1)
{
    A.obj->setkT(kT);
    stringstream ss;
    ss << kT;

    string base = "_kT=";
    base += ss.str();

    A.run_singlebond(1000000, 1000);
}

</code></pre>

Tetrahedral liquid:

<pre><code>

int n = 2000;
double packing_fraction = 0.01;

double l = cbrt(pi*(double)n/(6.*packing_fraction));

Condensate A(l, n);


TetrahedralPatch c(1.0, 1.4, 0.927);

//TwoTetrahedral c(10.0, 1.4, pi / 4., 0.0, 1., pi / 6., 0.0, 1., pi / 6., 1000, 1000);



BindingModelSingle b(0.99,0.01);

A.setBindingModel(b);

A.setpots(c);

//int a = system("python3 /home/dino/Documents/Condensate/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

//int a = system("python3 /home/dino/Desktop/tylercollab/Repo/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

//A.run_singlebond(10000, 1000);

A.setviscosity(1.0);


double beta = 1.0;

A.obj->setkT(1./beta);



stringstream ss;
ss << beta;

string base = "_beta=";
base += ss.str();

A.run_singlebond(10000, 1000, base);


</code></pre>

# Two Tetrahedrals and A Singleton Simulation

<pre><code>

int m1 = 250;
int m2 = m1+100;
BindingModelTernary b(m1*4,4*m2);


// b.setup(0.99,0.01,0.01,0.01,0.0,0.0,
// 0.,
// 0.0,
// 0.0,
// 0.0,
// 0.0,
// 0.0,
// 0.0,
// 0.0,
// 0.0);


b.setup(0.99,0.01,0.01,0.01,0.99,0.2,
0.,
8.,
0.0,
0.0,
0.0,
1.0,
0.0,
0.0,
0.0);


matrix<double> params(6,3);
params(0, 0) = 10.0; // 1<->1
params(0, 1) = 1.4;
params(0, 2) = 0.927;

params(1, 0) = 0.0; //1 <->2
params(1, 1) = 1.0;
params(1, 2) = 0.927;

params(2, 0) = 10.0; //1 <-> 3
params(2, 1) = 1.4;
params(2, 2) = 0.927;

params(3, 0) = 0.0; // 2 <-> 2
params(3, 1) = 1.0;
params(3, 2) = 0.927;

params(4, 0) = 10.0; // 2<-> 3
params(4, 1) = 1.4;
params(4, 2) = 0.927;

params(5, 0) = 0.0; //3 <-> 3
params(5, 1) = 1.0;
params(5, 2) = 0.927;







int n = 600;

TwoTetrahedralAndSingle c(params, m1, m2, n);
int n2 = 100;
double packing_fraction = 0.01;

double l = cbrt(pi*(double)n/(6.*packing_fraction));

Condensate A(l, n);


//TwoTetrahedral c(10.0, 1.4, pi / 4., 0.0, 1., pi / 6., 0.0, 1., pi / 6., 1000, 1000);

// string filp = "/home/dino/Desktop/Chemistry/SimulationResults/ChemicalOscillator/sim-20-12-14-19:43:58/pos_beta=1_i=0455.csv";
// string filo = "/home/dino/Desktop/Chemistry/SimulationResults/ChemicalOscillator/sim-20-12-14-19:43:58/orientation_beta=1_i=0455.csv";

// string filp = "/home/dino/Documents/Condensate/TernaryFluid2/pos_beta=1_i=02097.csv";
// string filo = "/home/dino/Documents/Condensate/TernaryFluid2/orientation_beta=1_i=02097.csv";

// double T;
// bool err1;
// bool err2;
// matrix<double> temppos = importcsv(filp, T, err1);
// matrix<double> tempori = importcsv(filo, T, err1);

// matrix<double> newpos(2000,3);
// matrix<double> newori(2000,3);

// A.obj->setdat(temppos);
// A.obj->setorientation(tempori);


// cout << b.tripr111 << endl;
// cout << b.doubr11 << endl;
// cout << b.doubr22 << endl;
// cout << b2.on_rate << endl;
// cout << b2.off_rate << endl;
// pausel();

//TetrahedralPatch c2(10.0,1.4,0.927);

A.setBindingModel(b);

A.setpots(c);

//int a = system("python3 /home/dino/Documents/Condensate/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

//int a = system("python3 /home/dino/Desktop/tylercollab/Repo/Code/Plotting/FigureMonitor.py ./ ./col.csv >filecreationlog &");

//A.run_singlebond(10000, 1000);

A.setviscosity(1.0);


double beta = 1.0;

A.obj->setkT(1./beta);



stringstream ss;
ss << beta;

string base = "_beta=";
base += ss.str();

A.run_singlebond(10000000, 1000, base);


</code></pre>

Analysis of clusters

<pre><code>

        #pragma omp parallel 
        {
            vector<mdpair> mypairs_private;
            mypairs_private.reserve(number_to_reserve);

            #pragma omp for nowait schedule(dynamic) 
            for (int i = 0; i < nbins.getsize() - 1; i++)
            {
                
                int size_of_cluster = nbins[i + 1] - nbins[i];
                //cout << "nbins: " << i << " " << size_of_cluster << endl;
                if (size_of_cluster == 1)
                {
                    //do nothing
                    int i1 = indexes[nbins[i]];
                    //isbound[i1]=false;

                    bo.isbound[i1] = false;

                    //not bound to anything.
                }
                else if (size_of_cluster == 2)
                {
                    //all fine, bindings
                    int ti1 = indexes[nbins[i]];
                    int ti2 = indexes[nbins[i] + 1];

                    int i1;
                    int i2;
                    sort_doublet(ti1, ti2, i1, i2);

                    bool alreadybound_to_eachother = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];

                    bool aft;

                    bm.doublet(alreadybound_to_eachother, i1, i2, aft);

                    if (aft)
                    {
                        bo.boundto[i1] = i2;
                        bo.boundto[i2] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = true;
                        //tbt++;
                        mypairs_private.push_back(mdpair(i1,i2));
                    }
                    else
                    {

                        bo.isbound[i1] = false;
                        bo.isbound[i2] = false;

                    }
                    //bool already_bound = prebound(i1,i2);
                }
                else if (size_of_cluster == 3)
                {

                    int ti1 = indexes[nbins[i]];
                    int ti2 = indexes[nbins[i] + 1];
                    int ti3 = indexes[nbins[i] + 2];

                    //SORT THE INDICES (IMPORTANT)

                    int i1;
                    int i2;
                    int i3;

                    sort_triplet(ti1, ti2, ti3, i1, i2, i3);


                    //DETERMINE WHETHER THEY ARE BOUND
                    bool b12 = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];
                    bool b23 = bo.boundto[i2] == i3 && bo.boundto[i3] == i2 && bo.isbound[i2] && bo.isbound[i3];
                    bool b13 = bo.boundto[i1] == i3 && bo.boundto[i3] == i1 && bo.isbound[i1] && bo.isbound[i3];

                    //DETERMINE THE CONNECTIVENESS OF THE GRAPH
                    //remember, that in order to count as a triplet
                    bool c12 = false;
                    bool c23 = false;
                    bool c13 = false;

                    int nb1 = tempbound[i1];

                    int nb2 = tempbound[i2];
                    int nb3 = tempbound[i3];

                    // if(nb1 > 2 || nb2 > 2 || nb3 > 2) {
                        
                    //     cout << i1 << " " << i2 << " " << i3 << endl;

                    //     outfunc(tempbound, "t1");
                    //     outfunc(boindices, "t2");

                    //     pausel();
                    //     error("something weird in code");
                    // }

                    if (nb1 == 1)
                    {
                        int tempi = boindices(i1, 0);
                        if (tempi == i2)
                        {
                            c12 = true;
                            c13 = false;
                            c23 = true; //in order to be a triplet
                        }
                        else if (tempi == i3)
                        {
                            c13 = true;
                            c12 = false;
                            c23 = true;
                        }
                        else
                            error("something weird");

                        //check the other
                    }
                    else if (nb1 == 2)
                    {
                        c12 = true;
                        c13 = true;

                        int nb2 = tempbound[i2];
                        if (nb1 == 1)
                        {
                            c23 = false;
                        }
                        else
                        {
                            c23 = true;
                        }
                    }
                    else
                    {
                        cout << size_of_cluster << endl;
                        cout << b12 <<  " " << b23 << " " << b13 << endl;
                        cout << i1 << " " << i2 << " " << i3 << endl;
                        cout << tempbound[i1] << " " << tempbound[i2] << " " << tempbound[i3] << endl;
                        
                        for(int k = 0 ; k < tempbound[i1] ; k++) {
                            cout << boindices(i1,k) << " ";
                        }
                        cout << endl;

                        for (int k = 0; k < tempbound[i2]; k++)
                        {
                            cout << boindices(i2, k) << " ";
                        }
                        cout << endl;

                        for (int k = 0; k < tempbound[i3]; k++)
                        {
                            cout << boindices(i3, k) << " ";
                        }
                        cout << endl;

                        error("error in clustering algorithm");
                    }

                    bool a12;
                    bool a23;
                    bool a13;
                    bm.triplet(b12, b23, b13, c12, c23, c13, i1, i2, i3, a12, a23, a13);

                    if (a12)
                    {
                        bo.boundto[i1] = i2;
                        bo.boundto[i2] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = true;
                        bo.isbound[i3] = false;

                        mypairs_private.push_back(mdpair(i2, i1));
                    }
                    else if (a23)
                    {
                        bo.boundto[i2] = i3;
                        bo.boundto[i3] = i2;
                        bo.isbound[i1] = false;
                        bo.isbound[i2] = true;
                        bo.isbound[i3] = true;

                        mypairs_private.push_back(mdpair(i3, i2));
                    }
                    else if (a13)
                    {
                        bo.boundto[i1] = i3;
                        bo.boundto[i3] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = false;
                        bo.isbound[i3] = true;

                        mypairs_private.push_back(mdpair(i3, i1));
                    }
                    else
                    {
                        bo.isbound[i1] = false;
                        bo.isbound[i2] = false;
                        bo.isbound[i3] = false;
                    }
                }
                else
                {
                    // cout << size_of_cluster << endl;
                    // for(int j = 0 ; j < size_of_cluster ; j++) {
                    //     for(int k = 0  ; k <tempbound[indexes[j]] ; k++ ) {
                    //         vector1<double> un(dimension);
                    //         double dis;
                    //         int p1,p2;
                    //         iny.which_particle(indexes[j], boindices(indexes[j],k), p1, p2);
                    //         geo->distance_vector(*dat, p1,p2, un, dis);


                    //         //un = i-j
                    //         dis = sqrt(dis);
                    //         cout << "{{" << indexes[j] << "\\" <<"[UndirectedEdge]" << boindices(indexes[j], k) <<"}," << dis << "},";
                                             
                    //     }
                    // }
                    // cout << endl;
                    // //for large clusters, repeat the partitioning 
                    // pausel();
                    if(size_of_cluster < 10) {
                    //firstly, obtain the graph of edges;
                    vector<mdpair> matched;
                    matched.reserve(size_of_cluster*(size_of_cluster-1)/2);
                    
                    for (int j = nbins[i]; j < nbins[i + 1]; j++)
                    {
                        for (int k = 0; k < tempbound[indexes[j]]; k++)
                        {
                            mdpair m1;
                            int f1 = indexes[j];
                            int f2 = boindices(indexes[j], k);
                            
                            if(f1<f2) {
                                m1.a = f1;
                                m1.b = f2;
                            }
                            else{
                                m1.a = f2;
                                m1.b = f1;
                            }
                            matched.push_back(m1);
                        }
                    }

                    

                    //next, remove duplicate edges

                    sort(matched.begin(), matched.end());
                    matched.erase(unique(matched.begin(), matched.end()), matched.end());

                    //cout << "removed duplicates" << endl;

                    //get the initial state
                    vector1<bool> bounded(matched.size());

                    for(int j = 0 ; j < matched.size() ; j++) {
                        int i1  = matched[j].a;
                        int i2  = matched[j].b;
                        bool b12 = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];
                        bounded[j] = b12;
                    }

                    //find all complimentary patterns

                    //cout << bounded << endl;
                    int m = 1 << matched.size();
                   // cout << m << endl;

                    if(abs(m)<1000) { //if the states are too big, don't do it

                    //get the set of all possible bools to move to
                    vector< vector1<bool> > possible_states;
                    
                    possible_states.reserve(m);

                    //cout << pow(2,matched.size()) << endl;
   

                    for(int j = 0 ; j < m ; j++) {
                        vector1<bool> binaryNum(matched.size());
                        int n = j;
                        int k = 0;
                        while(n > 0) {
                            binaryNum[k] = n % 2;
                            n = n/2;
                            k++;
                        }

                        if(IndependentEdge(matched,binaryNum)) {
                           // cout << binaryNum << endl;
                            possible_states.push_back(binaryNum);
                        }
                    }
                    
                    vector1<bool> afters(matched.size());
                    
                    bm.nlet(bounded, matched, possible_states, afters);


                    for(int j =  0 ; j < afters.getsize() ; j++) {
                        if(!afters[j]) {
                            int i1 = matched[j].a;
                            int i2 = matched[j].b;

                            bo.isbound[i1] = false;
                            bo.isbound[i2] = false;
                        }
                    }
                    // need to do it this way because of collisions

                    for (int j = 0; j < afters.getsize(); j++)
                    {
                        if (afters[j])
                        {
                            int i1 = matched[j].a;
                            int i2 = matched[j].b;

                            bo.isbound[i1] = true;
                            bo.isbound[i2] = true;
                            bo.boundto[i1] = i2;
                            bo.boundto[i2] = i1;

                            mypairs_private.push_back(mdpair(i1, i2));
                        }
                    }
                }
                    //cout << possible_states.capacity() << endl;
                    //cout << i1 << " " << size_of_cluster << endl;

                    //i1 is the number bound already
                // cout << bounded << endl;
                // for(int j  = 0  ; j < matched.size() ; j++) {
                // cout << matched[j] << endl;
                // }
                // cout << afters << endl;
                // pausel();
               // cout << "ncluster done" << endl;
                }
                }

            }

            #pragma omp for schedule(static) ordered
            for (int i = 0; i < omp_get_num_threads(); i++)
            {
            #pragma omp ordered
                mypairs.insert(mypairs.end(), mypairs_private.begin(), mypairs_private.end());
            }
        }

        </code></pre>


        Old large cluster algorithm

        <pre><code>
        {

                    if(size_of_cluster < 10) {
                    //firstly, obtain the graph of edges;
                    vector<mdpair> matched;
                    matched.reserve(size_of_cluster*(size_of_cluster-1)/2);
                    
                    for (int j = nbins[i]; j < nbins[i + 1]; j++)
                    {
                        //for (int k = 0; k < tempbound[indexes[j]]; k++)
                        for (int k = 0; k < tempbound[jhg[j].b]; k++)
                        {
                            mdpair m1;
                            //int f1 = indexes[j];
                            //int f2 = boindices(indexes[j], k);

                            int f1 = jhg[j].b;
                            int f2 = boindices(jhg[j].b, k);

                            if(f1<f2) {
                                m1.a = f1;
                                m1.b = f2;
                            }
                            else{
                                m1.a = f2;
                                m1.b = f1;
                            }
                            matched.push_back(m1);
                        }
                    }

                    //push back all possible pairs

                    

                    //next, remove duplicate edges

                    sort(matched.begin(), matched.end());
                    matched.erase(unique(matched.begin(), matched.end()), matched.end());

                    //cout << "removed duplicates" << endl;

                    //get the initial state
                    vector1<bool> bounded(matched.size());

                    for(int j = 0 ; j < matched.size() ; j++) {
                        int i1  = matched[j].a;
                        int i2  = matched[j].b;
                        bool b12 = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];
                        bounded[j] = b12;
                    }

                    //find all complimentary patterns

                    //cout << bounded << endl;
                    int m = 1 << matched.size();
                   // cout << m << endl;

                    if(abs(m)<1000) { //if the states are too big, don't do it

                    //get the set of all possible bools to move to
                    vector< vector1<bool> > possible_states;
                    
                    possible_states.reserve(m);

                    //cout << pow(2,matched.size()) << endl;
   

                    for(int j = 0 ; j < m ; j++) {
                        vector1<bool> binaryNum(matched.size());
                        int n = j;
                        int k = 0;
                        while(n > 0) {
                            binaryNum[k] = n % 2;
                            n = n/2;
                            k++;
                        }

                        if(IndependentEdge(matched,binaryNum)) {
                           // cout << binaryNum << endl;
                            possible_states.push_back(binaryNum);
                        }
                    }
                    
                    vector1<bool> afters(matched.size());
                    
                    bm.nlet(bounded, matched, possible_states, afters);


                    for(int j =  0 ; j < afters.getsize() ; j++) {
                        if(!afters[j]) {
                            int i1 = matched[j].a;
                            int i2 = matched[j].b;

                            bo.isbound[i1] = false;
                            bo.isbound[i2] = false;
                        }
                    }
                    // need to do it this way because of collisions

                    for (int j = 0; j < afters.getsize(); j++)
                    {
                        if (afters[j])
                        {
                            int i1 = matched[j].a;
                            int i2 = matched[j].b;

                            bo.isbound[i1] = true;
                            bo.isbound[i2] = true;
                            bo.boundto[i1] = i2;
                            bo.boundto[i2] = i1;

                            mypairs_private.push_back(mdpair(i1, i2));
                        }
                    }
                }
                    //cout << possible_states.capacity() << endl;
                    //cout << i1 << " " << size_of_cluster << endl;

                    //i1 is the number bound already
                // cout << bounded << endl;
                // for(int j  = 0  ; j < matched.size() ; j++) {
                // cout << matched[j] << endl;
                // }
                // cout << afters << endl;
                // pausel();
               // cout << "ncluster done" << endl;
                }
                }
                
        </code></pre>

Working parallel connected components

<pre><code>

void ConnectedComponentsParallel(matrix<int> &adj, vector1<int> &indexes) {
    vector1<int> f(indexes);
    vector1<int> fnext(indexes);
    //if(res.size() != indexes.size) error("error in capacity of save function");
    vector1<int> ftemp(indexes.getsize());

    // #pragma omp parallel for schedule(static)
    // for(int i = 0  ; i < indexes.getNsafe() ; i++) {
    //     f[i] =  indexes[i];
    //     gf[i] = indexes[i];
    // }
 

    int nr = adj.nrows;

    for(;;) {
        bool equiv;

        #pragma omp parallel for schedule(static)
        for(int i= 0 ; i < f.size ; i++) {
            ftemp.data[i] = f.data[i];
        }
        //#pragma omp parallel 
        
           // cout << omp_get_num_threads() << endl;
        #pragma omp parallel for schedule(static)
        for(int i = 0  ; i < nr ; i++) {
            int u = adj.mat[i*2 + 0];
            int v = adj.mat[i*2 + 1];
            // int u,v;
            // sort_doublet(u1,v1,u,v);
            
            int fu = f.data[u];
            //cout << u << " " << v << " " << fu << " " << f.data[v] << endl;
            if (fu == f.data[fu] && f.data[v] < fu)
            {
                fnext.data[fu] = f.data[v];
            }
            //pausel();
        }

        //cout << "loop done" << endl;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < f.size; i++)
            f.data[i] = fnext.data[i];         

        #pragma omp parallel for schedule(static)
        for(int i = 0 ; i < f.size ; i++) {
            int u = i;
            int fu = f.data[u];

            if (fu != f.data[fu])
            {
                fnext.data[u] = f.data[fu];
            }
        }

        equiv = true;
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < f.size; i++)
        {
            if (f.data[i] != fnext.data[i])
            {
                f.data[i] = fnext.data[i];
            }
            if(fnext.data[i] != ftemp.data[i]) {
                equiv = false;
            }
        }
        

        //cout << "done 3" << endl;    
        //bool equiv = true;
       if(equiv) break;


    }
    
    indexes = f;
    // cout << f << endl;
    // pausel();
    // //f = gf;
    // //res.resize(indexes.getsize());
    // //vector<mdpair> res(indexes.getsize());
    // #pragma omp parallel for
    // for(int i = 0 ; i < f.getsize() ; i++) {
    //     res[i].a = f.data[i];
    //     res[i].b = indexes.data[i];
    // }
    // mdpairsort fg;
    // std::sort(res.begin(),res.end(),fg);
    //return res;

}
</code></pre>

Bivalent run

<pre><code>
srand (time(NULL));
double packing_fraction;
double int1;
double int2;
double int3;
int runtime;
if(argc==6) {
runtime = atof(argv[1]);
packing_fraction = atof(argv[2]);
int1 = atof(argv[3]);
int2 = atof(argv[4]);
int3 = atof(argv[5]);
}
else {
    runtime = 10000;
    packing_fraction = 0.01;
    int1 =10.0;
    int2 = 10.0;
    int3 = 10.0;
}
matrix<double> params(28,3);

cout << packing_fraction << " " << int1 << " " << int2 << " " << int3 << endl;

// for(int i = 0 ; i < 28 ; i++) {
//     params(i, 0) =  10.0; //strength 
//     params(i, 1) =  1.4; //distance
//     params(i, 2) = 0.927;
// }

int iter = 0;
for(int i = 0  ; i < 4 ; i++) {
    for(int j = 0  ; j < 4 ; j++) {
        params(iter,0) = int1;
        params(iter,1) = 1.4;
        params(iter,2) =  0.927;
        iter++;
    }
}

for(int i = 0  ; i < 4 ; i++) {
    for(int j = 0 ; j < 2 ; j++) {
        if(j==1) {
            params(iter, 0) = int2;
            params(iter, 1) = 1.4*0.75;
            params(iter, 2) = 1.1;
        }
        else{
            params(iter, 0) = 0.0;
            params(iter, 1) = 1.0 * 0.75;
            params(iter, 2) = 0.927;
        }
        iter++;
    }
}

for (int i = 0; i < 2; i++)
{
    for (int j = 0; j < 2; j++)
    {

        if(i==0 && j ==0) {
            params(iter, 0) = int3;
            params(iter, 1) = 1.4*0.5;
            params(iter, 2) = 1.1;
        }
        else{
            params(iter, 0) = 0.0;
            params(iter, 1) = 1. * 0.5;
            params(iter, 2) = 0.927;
        }

        iter++;
    }
}


int n = 10000;
int nt = 2000;
TetrahedralWithBivalent c(params,nt,n);

TetrahedralWithSingle c2(10.0, 1.4, 0.927, 10., 1.4, 0.927, 10.0, 1.4, 0.927, nt, n);


TetrahedralPatch c3(10.0,1.4,0.927);

//double packing_fraction = 0.02;

double l = cbrt(pi * (double)n / (6. * packing_fraction));

BindingModelBinary b(nt*4);

b.setup_equilibrium();

Condensate A(l, n);

A.setBindingModel(b);

A.setpots(c);

A.setviscosity(0.1);

double beta = 1.;

A.obj->setkT(1. / beta);

stringstream ss;
ss << beta;

string base = "_beta=";
base += ss.str();

A.run_singlebond_different_sizes(runtime, 1000, nt, base);
</code></pre>

Langevin Force

<pre><code>

void LangevinNVTR::calculate_forces_and_torques3D_onlyone(matrix<int> &pairs, ComboPatch &iny, BinaryBindStore &bo, AbstractBindingModel &bm, matrix<double> &forces, matrix<double> &torques)
{

    //for all the pairs, for all bindings

    //for a given sphere geometry

    //int np1 = sqrt(iny.getsize());
    int total_number_of_patches = bo.boundto.getsize();//iny.get_total_patches(this->getN());


    vector1<int> tempbound(total_number_of_patches,0); //no binding to begin wtih

    int depth_of_matrix = 10; //Choose this value to be deep enough such that all values can be stored

    matrix<int> boindices(total_number_of_patches, depth_of_matrix);

    // vector<mdpair> edgelist;
    // edgelist.reserve(total_number_of_patches);

    std::mutex mtx;

    unsigned int i;

    #pragma omp parallel for schedule(dynamic)
    for (i = 0; i < pairs.getNsafe(); ++i)
    {
        int p1 = pairs(i, 0);
        int p2 = pairs(i, 1);
        if(p2 < p1) { //INDICES NEED TO BE SORTED FOR IT TO WORK
            int tp1 = p1;
            p1 = p2;
            p2 = tp1;

        }
        //int i1 = pairs(i,2);
        double dis;
        //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
        vector1<double> un(dimension);
        geo->distance_vector(*dat, p1, p2, un, dis);

        //un = i-j
        dis = sqrt(dis);

        un /= dis;
        double dx = un.gpcons(0);
        double dy = un.gpcons(1);
        double dz = un.gpcons(2);

        double qtemp0 = orient->gpcons(p1, 0);
        double qtemp1 = orient->gpcons(p1, 1);
        double qtemp2 = orient->gpcons(p1, 2);
        double qtemp3 = orient->gpcons(p1, 3);
        double qtemp4 = orient->gpcons(p1, 4);
        double qtemp5 = orient->gpcons(p1, 5);
        double qtemp6 = orient->gpcons(p1, 6);
        double qtemp7 = orient->gpcons(p1, 7);
        double qtemp8 = orient->gpcons(p1, 8);

        double gtemp0 = orient->gpcons(p2, 0);
        double gtemp1 = orient->gpcons(p2, 1);
        double gtemp2 = orient->gpcons(p2, 2);
        double gtemp3 = orient->gpcons(p2, 3);
        double gtemp4 = orient->gpcons(p2, 4);
        double gtemp5 = orient->gpcons(p2, 5);
        double gtemp6 = orient->gpcons(p2, 6);
        double gtemp7 = orient->gpcons(p2, 7);
        double gtemp8 = orient->gpcons(p2, 8);

        // for (int j = 0; j < iny.num_patches(p1) ; j++)
        // {
        //     for (int k = 0; k < iny.num_patches(p2); k++)
        //     {

                //int potn = np1 * j + k;
                
            int **q = new int*;
            if(iny.safe) {
            iny.UpdateIterator(p1,p2);
            *q= *iny.p;
            }
            else{
            iny.UpdateIteratorSafe(p1,p2,q);

            }


            //int **q = iny.p;

            for (int tp = 1; tp < (*q)[0] + 1; tp++)
            {
                int potn = (*q)[tp];
                // vector1<double> params = (iny.potential_bundle)[potn]->getparameters();
                //cout << potn << endl;
                double nxb1;// = params[0]; //iny[potn]->nxb1;
                double nyb1;// = params[1]; //iny[potn]->nyb1;
                double nzb1;// = params[2]; //iny[potn]->nzb1;

                double nxb2;// = params[3]; //iny[potn]->nxb2;
                double nyb2;// = params[4]; //iny[potn]->nyb2;
                double nzb2;// = params[5]; //iny[potn]->nzb2;

                double disp;// = params[6];

                double thetam;// = params[8];

//                cout << p1 << " " << p2 << " " << potn << endl;
                iny.get_params(p1,p2,potn,nxb1,nyb1,nzb1,nxb2,nyb2,nzb2,disp,thetam); //for this potential, get all the parameters

                double nx1 = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
                double ny1 = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
                double nz1 = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;

                double nx2 = nxb2 * gtemp0 + nyb2 * gtemp3 + nzb2 * gtemp6;
                double ny2 = nxb2 * gtemp1 + nyb2 * gtemp4 + nzb2 * gtemp7;
                double nz2 = nxb2 * gtemp2 + nyb2 * gtemp5 + nzb2 * gtemp8;

                double argthetai = -(nx1 * dx + ny1 * dy + nz1 * dz);
                double argthetaj = (nx2 * dx + ny2 * dy + nz2 * dz);

                int wp1, wp2;
                iny.which_patch(p1, p2, potn, wp1, wp2);

                // cout << p1 << " " << p2 << " " << wp1 << " " << wp2 << " " << disp << " " << thetam << endl;

                //different conditions depending on whether there is binding or not.
                double disp2;
                bool cond1 =  bo.boundto[wp1] == wp2;
                bool b1,b2;
                b1 = bo.isbound[wp1];
                b2 = bo.isbound[wp2];

                if(b1 && b2 && cond1) { //both bound and to each other
                    disp2 = disp;
                }
                else if (b1  && b2 && !cond1) //both bound and not to each other
                {
                    disp2 = 0.8*disp; //if both bound, make the conditions more onerous
                }
                else if (!b1 != !b2) //only one bound
                {
                    disp2 = 0.8*disp; //more onerous
                }
                else{
                    //neither bound
                    disp2 = disp;

                }


                if (argthetai > cos(thetam) && argthetaj > cos(thetam) && dis < disp2)
                {
                    //cout << argthetai <<  " " << argthetaj << endl;
                    // cout << p1 << " " << p2 << endl;
                    // cout << "params: " << endl;
                    // cout << qtemp0 << " " << qtemp1 << " " << qtemp2 << " " << qtemp3 << " " << qtemp4 << " " << qtemp5 << " " << qtemp6 << " " << qtemp7 << " " << qtemp8 << endl;
                    // cout << gtemp0 << " " << gtemp1 << " " << gtemp2 << " " << gtemp3 << " " << gtemp4 << " " << gtemp5 << " " << gtemp6 << " " << gtemp7 << " " << gtemp8 << endl;
                    // cout << un.gpcons(0) << " " << un.gpcons(1) << " " << un.gpcons(2) << endl;
                    // cout << nxb1 << " " << nyb1 << " " << nzb1 << endl;
                    // cout << nxb2 << " " << nyb2 << " " << nzb2 << endl;
                    // cout << nx1 << " " << ny1 << " " << nz1 << endl;
                    // cout << nx2 << " " << ny2 << " " << nz2 << endl;
                    // cout << argthetai << " " << argthetaj << " " << cos(thetam) << endl;
                    // cout << endl;
                    // if((wp1 > 400 &&wp2 > 400)) {
                    //     cout << "\n\n\n";
                    //     cout << "possible patch" << endl;
                    //     cout << p1 << " " << p2 << endl;
                    //     cout << wp1 << " " << wp2 << endl;
                    //     cout <<  dis << " " << disp << endl;
                    //     cout << argthetai << " " << argthetaj << " " <<cos(thetam) << endl;
                    //     cout << potn << endl;
                    //     cout << bo.isbound[wp1] << " " << bo.isbound[wp2] << endl;
                    //     cout << bo.boundto[wp1] << " " << bo.boundto[wp2] << endl;
                    //     cout << "\n\n\n next";
                    //     pausel();
                    // }
                   // cout << potn << endl;
                    
                    
                    //#pragma omp atomic read
                    //mtx.lock();

                    const std::lock_guard<std::mutex> lock(mtx);
                    
                    int iterator1 = tempbound[wp1];
                    //#pragma omp atomic read
                    int iterator2 = tempbound[wp2];


                    tempbound[wp1]++;
                    tempbound[wp2]++;

                    //mtx.unlock();
                   
                    boindices(wp1, iterator1) = wp2;
                    
                    boindices(wp2, iterator2) = wp1;


                }
            }
            //pausel();
            delete q;
            
        //     }
        // }
        }

        // cout << "which patch" << endl;
        //pausel();
    

    // for(int i =  0 ; i <total_number_of_patches ; i++) {
    //     for(int j = 0 ; j < tempbound[i] ; j++) {
    //         if(tempbound[i] >0 && tempbound[boindices(i,j)] == 0) {
                
    //             int bk = boindices(i,j);
    //             cout << "parts" << endl;
    //             cout << i << endl;
    //             cout << bk << endl;
    //             cout << "lens" << endl;
    //             cout << tempbound[i] << endl;
    //             cout << tempbound[bk] << endl;

    //             for(int l1 = 0  ; l1 < tempbound[i] ; l1++)
    //                 cout << boindices(i,l1) <<  " ";
    //             cout << endl;

    //             for(int l2 = 0  ; l2 < tempbound[bk] ; l2++)
    //                 cout << boindices(bk, l2) << " ";
    //             cout << endl;

    //             for (int l1 = 0; l1 < depth_of_matrix; l1++)
    //                 cout << boindices(i, l1) << " ";
    //             cout << endl;

    //             for (int l2 = 0; l2 < depth_of_matrix; l2++)
    //                 cout << boindices(bk, l2) << " ";
    //             cout << endl;

    //             error("asymmetry in bond list");
    //         }
    //     }
    // }

// //   

    #pragma omp parallel for schedule(static)
    for(int i = 0  ; i <total_number_of_patches ; i++) { //check bindings, if distance metric is wrong, break bindings
        if(bo.isbound[i]) { //if it is bound
            int bt = bo.boundto[i]; //it is bound to what

            int outside = true;
            for(int k = 0 ; k < tempbound[i] ; k++) {
                int pbt  = boindices(i,k);
                if(bt  == pbt ) {
                    outside = false;
                    break;
                }
            }
            if(outside) {
                bo.isbound[i]= false;
             //   cout << "unbind event" << endl;
               // cout << bo.isbound << endl;
               
                //cout << bo.isbound << endl;
            }
            
        }
    }


    matrix<int> edgelist = this->CreateEdgeList(boindices, tempbound);

    string sg = "a";
    vector1<int> indexes2(total_number_of_patches,sg);
    //std::vector<mdpair> jhg(total_number_of_patches);
    
    ConnectedComponentsParallel(edgelist, indexes2);
    
    matrix<int> boindices2(total_number_of_patches, depth_of_matrix);
    vector1<int> ccs(total_number_of_patches);

    
    // //#pragma omp parallel for schedule(static)
    for (int i = 0; i < total_number_of_patches; ++i)
    {
        int wp1 = indexes2[i];
        //mtx.lock();
        //const std::lock_guard<std::mutex> lock(mtx);

        int iterator1 = ccs[wp1];
        if(iterator1 < 10) {
            
            boindices2(wp1, iterator1) = i;
            ccs[wp1]++;
            // mtx.unlock();
        }
        //mtx.unlock();
    }
    
    

        int number_to_reserve = MIN(2*((total_number_of_patches+1)- total_number_of_patches),total_number_of_patches/2 );
//        // cout << number_to_reserve << endl;
        vector<mdpair> mypairs;//(number_to_reserve);
        mypairs.reserve(number_to_reserve);


        #pragma omp parallel 
        {
            vector<mdpair> mypairs_private;
            mypairs_private.reserve(number_to_reserve);

            #pragma omp for nowait schedule(dynamic) 
            for (int i = 0; i < total_number_of_patches; i++)
            {
                int size_of_cluster = ccs[i];

                if(size_of_cluster ==0) {

                }

                else if (size_of_cluster == 1)
                {
                    //do nothing
                    //int i1 = indexes[nbins[i]];
                    int i1 = boindices2(i,0);
                    //isbound[i1]=false;

                    bo.isbound[i1] = false;

                    //not bound to anything.
                }
                else if (size_of_cluster == 2)
                {
                    //all fine, bindings
                    // int ti1 = indexes[nbins[i]];
                    // int ti2 = indexes[nbins[i] + 1];

                    int ti1 = boindices2(i, 0);
                    int ti2 = boindices2(i, 1);

                    int i1;
                    int i2;
                    sort_doublet(ti1, ti2, i1, i2);

                    bool alreadybound_to_eachother = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];

                    bool aft;

                    bm.doublet(alreadybound_to_eachother, i1, i2, aft);

                    if (aft)
                    {
                        bo.boundto[i1] = i2;
                        bo.boundto[i2] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = true;
                        //tbt++;
                        mypairs_private.push_back(mdpair(i1,i2));
                    }
                    else
                    {

                        bo.isbound[i1] = false;
                        bo.isbound[i2] = false;

                    }
                    //bool already_bound = prebound(i1,i2);
                }
                else if (size_of_cluster == 3)
                {

                    // int ti1 = indexes[nbins[i]];
                    // int ti2 = indexes[nbins[i] + 1];
                    // int ti3 = indexes[nbins[i] + 2];

                    int ti1 = boindices2(i, 0);
                    int ti2 = boindices2(i, 1);
                    int ti3 = boindices2(i, 2);
                    //SORT THE INDICES (IMPORTANT)

                    int i1;
                    int i2;
                    int i3;

                    sort_triplet(ti1, ti2, ti3, i1, i2, i3);


                    //DETERMINE WHETHER THEY ARE BOUND
                    bool b12 = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];
                    bool b23 = bo.boundto[i2] == i3 && bo.boundto[i3] == i2 && bo.isbound[i2] && bo.isbound[i3];
                    bool b13 = bo.boundto[i1] == i3 && bo.boundto[i3] == i1 && bo.isbound[i1] && bo.isbound[i3];

                    //DETERMINE THE CONNECTIVENESS OF THE GRAPH
                    //remember, that in order to count as a triplet
                    bool c12 = false;
                    bool c23 = false;
                    bool c13 = false;

                    int nb1 = tempbound[i1];

                    int nb2 = tempbound[i2];
                    int nb3 = tempbound[i3];

                    // if(nb1 > 2 || nb2 > 2 || nb3 > 2) {
                        
                    //     cout << i1 << " " << i2 << " " << i3 << endl;

                    //     outfunc(tempbound, "t1");
                    //     outfunc(boindices, "t2");

                    //     pausel();
                    //     error("something weird in code");
                    // }

                    if (nb1 == 1)
                    {
                        int tempi = boindices(i1, 0);
                        if (tempi == i2)
                        {
                            c12 = true;
                            c13 = false;
                            c23 = true; //in order to be a triplet
                        }
                        else if (tempi == i3)
                        {
                            c13 = true;
                            c12 = false;
                            c23 = true;
                        }
                        else
                            error("something weird");

                        //check the other
                    }
                    else if (nb1 == 2)
                    {
                        c12 = true;
                        c13 = true;

                        int nb2 = tempbound[i2];
                        if (nb1 == 1)
                        {
                            c23 = false;
                        }
                        else
                        {
                            c23 = true;
                        }
                    }
                    else
                    {
                        cout << ccs[i] << " " << endl;
                        cout << size_of_cluster << endl;
                        cout << b12 <<  " " << b23 << " " << b13 << endl;
                        cout << i1 << " " << i2 << " " << i3 << endl;
                        cout << tempbound[i1] << " " << tempbound[i2] << " " << tempbound[i3] << endl;
                        
                        cout << indexes2[i1] << " " << indexes2[i2] << " " << indexes2[i3] << endl;

                        for(int k = 0  ; k < total_number_of_patches ; k++) {
                            if(indexes2[k] ==indexes2[i1]) {
                                cout << k << " ";
                            }
                        }
                        cout << endl;
                        cout << "indexes done" << endl;
                        cout << "parallel vs serial" << endl;
                        for(int k = 0 ; k < depth_of_matrix ; k++) {
                        cout << boindices2(i,k) <<  endl;
                        }
                        cout << endl;

                        for(int k = 0 ; k < depth_of_matrix ; k++) {
                            cout << boindices(i1, k) << " " << indexes2[boindices(i1, k)] << endl;
                            // cout << indexes2[boindices(i1, k)] << endl;
                        }
                        cout << endl;

                        for (int k = 0; k < depth_of_matrix; k++)
                        {
                            cout << boindices(i2, k) << " " << indexes2[boindices(i2, k)] << endl;
                        }
                        cout << endl;

                        for (int k = 0; k < depth_of_matrix; k++)
                        {
                            cout << boindices(i3, k) <<  " " << indexes2[boindices(i3, k)] << endl;
                        }


                        cout << endl;

                        error("error in clustering algorithm");
                    }

                    bool a12;
                    bool a23;
                    bool a13;
                    bm.triplet(b12, b23, b13, c12, c23, c13, i1, i2, i3, a12, a23, a13);

                    if (a12)
                    {
                        bo.boundto[i1] = i2;
                        bo.boundto[i2] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = true;
                        bo.isbound[i3] = false;

                        mypairs_private.push_back(mdpair(i2, i1));
                    }
                    else if (a23)
                    {
                        bo.boundto[i2] = i3;
                        bo.boundto[i3] = i2;
                        bo.isbound[i1] = false;
                        bo.isbound[i2] = true;
                        bo.isbound[i3] = true;

                        mypairs_private.push_back(mdpair(i3, i2));
                    }
                    else if (a13)
                    {
                        bo.boundto[i1] = i3;
                        bo.boundto[i3] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = false;
                        bo.isbound[i3] = true;

                        mypairs_private.push_back(mdpair(i3, i1));
                    }
                    else
                    {
                        bo.isbound[i1] = false;
                        bo.isbound[i2] = false;
                        bo.isbound[i3] = false;
                    }
                }
                else
                {

                    if(size_of_cluster < 10) {
                    //firstly, obtain the graph of edges;
                    vector<mdpair> matched;
                    matched.reserve(size_of_cluster*(size_of_cluster-1)/2);
                    
                    for (int j = 0; j < size_of_cluster; j++)
                    {
                        int myindex = boindices2(i, j);
                            //for (int k = 0; k < tempbound[indexes[j]]; k++)
                            for (int k = 0; k < tempbound[myindex]; k++)
                            {
                                mdpair m1;
                                //int f1 = indexes[j];
                                //int f2 = boindices(indexes[j], k);

                                int f1 = myindex;
                                int f2 = boindices(myindex, k);

                                if(f1<f2) {
                                    m1.a = f1;
                                    m1.b = f2;
                                }
                                else{
                                    m1.a = f2;
                                    m1.b = f1;
                                }
                                matched.push_back(m1);
                            }
                    }

                    //push back all possible pairs

                    

                    //next, remove duplicate edges

                    sort(matched.begin(), matched.end());
                    matched.erase(unique(matched.begin(), matched.end()), matched.end());

                    //cout << "removed duplicates" << endl;

                    //get the initial state
                    vector1<bool> bounded(matched.size());

                    for(int j = 0 ; j < matched.size() ; j++) {
                        int i1  = matched[j].a;
                        int i2  = matched[j].b;
                        bool b12 = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];
                        bounded[j] = b12;
                    }

                    //find all complimentary patterns

                    //cout << bounded << endl;
                    int m = 1 << matched.size();
                   // cout << m << endl;

                    if(abs(m)<1000) { //if the states are too big, don't do it

                    //get the set of all possible bools to move to
                    vector< vector1<bool> > possible_states;
                    
                    possible_states.reserve(m);

                    //cout << pow(2,matched.size()) << endl;
   

                    for(int j = 0 ; j < m ; j++) {
                        vector1<bool> binaryNum(matched.size());
                        int n = j;
                        int k = 0;
                        while(n > 0) {
                            binaryNum[k] = n % 2;
                            n = n/2;
                            k++;
                        }

                        if(IndependentEdge(matched,binaryNum)) {
                           // cout << binaryNum << endl;
                            possible_states.push_back(binaryNum);
                        }
                    }
                    
                    vector1<bool> afters(matched.size());
                    
                    bm.nlet(bounded, matched, possible_states, afters);


                    for(int j =  0 ; j < afters.getsize() ; j++) {
                        if(!afters[j]) {
                            int i1 = matched[j].a;
                            int i2 = matched[j].b;

                            bo.isbound[i1] = false;
                            bo.isbound[i2] = false;
                        }
                    }
                    // need to do it this way because of collisions

                    for (int j = 0; j < afters.getsize(); j++)
                    {
                        if (afters[j])
                        {
                            int i1 = matched[j].a;
                            int i2 = matched[j].b;

                            bo.isbound[i1] = true;
                            bo.isbound[i2] = true;
                            bo.boundto[i1] = i2;
                            bo.boundto[i2] = i1;

                            mypairs_private.push_back(mdpair(i1, i2));
                        }
                    }
                }
                    //cout << possible_states.capacity() << endl;
                    //cout << i1 << " " << size_of_cluster << endl;

                    //i1 is the number bound already
                // cout << bounded << endl;
                // for(int j  = 0  ; j < matched.size() ; j++) {
                // cout << matched[j] << endl;
                // }
                // cout << afters << endl;
                // pausel();
               // cout << "ncluster done" << endl;
                }
                }

            }

            #pragma omp for schedule(static) ordered
            for (int i = 0; i < omp_get_num_threads(); i++)
            {
            #pragma omp ordered
                mypairs.insert(mypairs.end(), mypairs_private.begin(), mypairs_private.end());
            }
        }

        

//     //    cout << "calc patches" << endl;

//         //Having gone through all the patches to determine the bindings, we can now calculate the forces

// /*         vector1<bool> visited(total_number_of_patches); //keep track of patches that we already calculated
//         vector1<int> par1(tbt);
//         vector1<int> par2(tbt);
//         int tb = 0;
//         int tb2 = 0;
//         int iter3 = 0;
//         for (int i = 0; i < total_number_of_patches; i++)
//         {
//             tb2 += bo.isbound[i];
//             if (bo.isbound[i] == true &&  visited[i] == false)
//             {
//                 visited[i] = true;
//                 visited[bo.boundto[i]] = true;
//                 //tb += bo.isbound[i];
//                 par1[iter3] = i;
//                 par2[iter3] = bo.boundto[i];
//                 iter3++;
//             }
//         }

//         cout << par1 << endl;
//         cout << par2 << endl;
//         for(int j  = 0  ; j < mypairs.size() ; j++) {
//             cout << mypairs[j].a << ", ";
//         }
//         cout << endl;
//         for (int j = 0; j < mypairs.size(); j++)
//         {
//             cout << mypairs[j].b << ", ";
//         }
//         cout << endl;

//         cout << mypairs.size() << endl;

//         cout << tbt << endl;
//         cout << tb2 << endl;
//         cout << tb << endl;
//         pausel(); */




        #pragma omp parallel for
        for(int i = 0  ; i < mypairs.size() ; i++) {
            int p1;
            int p2;
            int wp1 = mypairs[i].a;
            int wp2 = mypairs[i].b;
            iny.which_particle(wp1,wp2,p1,p2);
            
            if (p2 < p1)
            { //INDICES NEED TO BE SORTED FOR IT TO WORK
                int tp1 = p1;
                p1 = p2;
                p2 = tp1;
            }
            double dis;
            vector1<double> un(dimension);
            geo->distance_vector(*dat, p1, p2, un, dis);

            int potn = iny.which_potential(p1, p2, wp1, wp2);

            dis = sqrt(dis);

            un /= dis;




            double fx;
            double fy;
            double fz;

            double tix;
            double tiy;
            double tiz;

            double tjx;
            double tjy;
            double tjz;

           // cout << p1 << " " << p2 << endl;
            (iny.potential_bundle)[potn]->force_and_torque(un, dis, *orient, p1, p2, fx, fy, fz, tix, tiy, tiz, tjx, tjy, tjz);
;
            /* 
            double dx = un.gpcons(0);
            double dy = un.gpcons(1);
            double dz = un.gpcons(2);

            double qtemp0 = orient->gpcons(p1, 0);
            double qtemp1 = orient->gpcons(p1, 1);
            double qtemp2 = orient->gpcons(p1, 2);
            double qtemp3 = orient->gpcons(p1, 3);
            double qtemp4 = orient->gpcons(p1, 4);
            double qtemp5 = orient->gpcons(p1, 5);
            double qtemp6 = orient->gpcons(p1, 6);
            double qtemp7 = orient->gpcons(p1, 7);
            double qtemp8 = orient->gpcons(p1, 8);

            double gtemp0 = orient->gpcons(p2, 0);
            double gtemp1 = orient->gpcons(p2, 1);
            double gtemp2 = orient->gpcons(p2, 2);
            double gtemp3 = orient->gpcons(p2, 3);
            double gtemp4 = orient->gpcons(p2, 4);
            double gtemp5 = orient->gpcons(p2, 5);
            double gtemp6 = orient->gpcons(p2, 6);
            double gtemp7 = orient->gpcons(p2, 7);
            double gtemp8 = orient->gpcons(p2, 8);

           
            double nxb1; // = params[0]; //iny[potn]->nxb1;
            double nyb1; // = params[1]; //iny[potn]->nyb1;
            double nzb1; // = params[2]; //iny[potn]->nzb1;

            double nxb2; // = params[3]; //iny[potn]->nxb2;
            double nyb2; // = params[4]; //iny[potn]->nyb2;
            double nzb2; // = params[5]; //iny[potn]->nzb2;

            double disp; // = params[6];

            double thetam; // = params[8];

            //                cout << p1 << " " << p2 << " " << potn << endl;
            iny.get_params(p1, p2, potn, nxb1, nyb1, nzb1, nxb2, nyb2, nzb2, disp, thetam); //for this potential, get all the parameters

            double nx1 = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
            double ny1 = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
            double nz1 = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;

            double nx2 = nxb2 * gtemp0 + nyb2 * gtemp3 + nzb2 * gtemp6;
            double ny2 = nxb2 * gtemp1 + nyb2 * gtemp4 + nzb2 * gtemp7;
            double nz2 = nxb2 * gtemp2 + nyb2 * gtemp5 + nzb2 * gtemp8;

            double argthetai = -(nx1 * dx + ny1 * dy + nz1 * dz);
            double argthetaj = (nx2 * dx + ny2 * dy + nz2 * dz);
            cout << (iny.potential_bundle)[potn]->getparameters() << endl;
            cout << "particles:  " <<p1 << " " << p2 << endl;
            cout << "distance: " << dis << endl;
            cout << "potential number: " << potn << endl;
            cout << "which patch: " << wp1 << " " << wp2 << endl;
            cout << "force: " << fx << " " << fy << " " << fz << endl;
            cout << "torques: " << tix << " " << tiy << " " << tiz << endl;
            cout << "orientation1: " << qtemp0 << " " << qtemp1 << " " << qtemp2 << " " << qtemp3 << " " << qtemp4 << " " << qtemp5 << " " << qtemp6 << " " << qtemp7 << " " << qtemp8 << endl;
            cout << "orientation2: " << gtemp0 << " " << gtemp1 << " " << gtemp2 << " " << gtemp3 << " " << gtemp4 << " " << gtemp5 << " " << gtemp6 << " " << gtemp7 << " " << gtemp8 << endl;
            cout << "dx: " << dx << " " << dy << " " << dz << endl;
            cout << "n: " << nx1 << " " << ny1 << " " << nz1 << endl;
            cout << "n2: " << nx2 << " " << ny2 << " " << nz2 << endl;
            cout << "arguments: " << argthetai << " " << argthetaj << " " << cos(thetam) << endl;
            pausel();  */

            forces(p1, 0) += fx;
            forces(p1, 1) += fy;
            forces(p1, 2) += fz;

            forces(p2, 0) += -fx;
            forces(p2, 1) += -fy;
            forces(p2, 2) += -fz;

            torques(p1, 0) += tix;
            torques(p1, 1) += tiy;
            torques(p1, 2) += tiz;

            torques(p2, 0) += tjx; // - dis * (fz * un[1] - fy * un[2]);
            torques(p2, 1) += tjy; // - dis * (fz * un[0] + fx * un[2]);
            torques(p2, 2) += tjz; // - dis * (fy * un[0] - fx * un[1]);
        }
//UP TO HERE

        //check std atomic
        
       // cout << "forces" << endl;
        /* 
        #pragma omp parallel for
        for (int i = 0; i < total_number_of_patches; ++i)
        {
            if (bo.isbound[i] == true && visited[i] == false) //only for bound patches we haven't visisted do we calculate forces
            {

                visited[i] = true;             // We have now visisted this patch
                visited[bo.boundto[i]] = true; //we have also visisted the patch that it is bound to

                // int p1 = floor(i / np1);             //particle number 1
                // int p2 = floor(bo.boundto[i] / np1); //particle number 2

                //get the particle numbers from the patch numbers:

                int p1;
                int p2;

                iny.which_particle(i, bo.boundto[i], p1, p2);

                double dis;
                //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
                vector1<double> un(dimension);
                geo->distance_vector(*dat, p1, p2, un, dis);

                //un = i-j

                int potn = iny.which_potential(p1, p2, i, bo.boundto[i]);
                // int potn = (i % np1) * np1 + (bo.boundto[i] % np1);
                dis = sqrt(dis);

                un /= dis;

                double fx;
                double fy;
                double fz;

                double tix;
                double tiy;
                double tiz;

                double tjx;
                double tjy;
                double tjz;

                (iny.potential_bundle)[potn]->force_and_torque(un, dis, *orient, p1, p2, fx, fy, fz, tix, tiy, tiz, tjx, tjy, tjz);


                forces(p1, 0) += fx;
                forces(p1, 1) += fy;
                forces(p1, 2) += fz;

                forces(p2, 0) += -fx;
                forces(p2, 1) += -fy;
                forces(p2, 2) += -fz;

                torques(p1, 0) += tix;
                torques(p1, 1) += tiy;
                torques(p1, 2) += tiz;

                torques(p2, 0) += tjx; // - dis * (fz * un[1] - fy * un[2]);
                torques(p2, 1) += tjy; // - dis * (fz * un[0] + fx * un[2]);
                torques(p2, 2) += tjz; // - dis * (fy * un[0] - fx * un[1]);
            }
            else
            {
                //do nothing
            }
        } */
   // cout << "calc forces" << endl;

    //now we have only the real forces, we no longer need to calculate the forces for the non-bound particles:
    }
    
</code></pre>