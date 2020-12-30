In this directory, there exist various classes that define basic data types necessary for mathematical applications, such as vectors and matrices. The header files vector1.h and matrix2.h define the methods that operate on these classes. In this document, some simple examples of the operations of these classes will be given.

Both of these classes use templates, and thus we can call a vector or a matrix of any class in which we are interested. For example, to initialize a vector of length 3, where the data type is doubles, we would write:

```
vector1<double> a(3)
```
This would initialize a vector called a filled with doubles. The default setting is to initialize all the values of the vector as 0. I.e., running the following

```
vector1<double> a(3)
cout << a;
```
would lead to an output
```ruby
{0.,0.,0.}
```

If we wish to initialize with a value that is not zero (e.g 1.) we would instead call:

```
vector1<double> a(3,1.)
```

Once we have created the vector object, we can access its elements using square brackets, for example in this code snipped we create a vector, output one of the values, then change that value and output again:

```
vector1<double> a(3);
cout << a[0] << "\n";
a[0] = 1.;
cout << a[0] << "\n";
```

which would lead to the output:
```
0.
1.
```

Matrices can be defined in a similar way, except now we have an additional parameter accounting for both the rows and the columns of the matrix, for example creating a 3 by 3 matrix of doubles would look like:
```
matrix<double> a(3,3)
```

where once again each element of the matrix is initialized to zero. If we wish to access or modify elements of matrices, we use the curly bracket operator, for example in the following we initialize a matrix, and change the (0,0) index to be a different value:

```
matrix<double> a(3,3);
cout << a(0,0) << "\n";
a(0,0) = 1.;
cout << a(0,0) << "\n";
```

where again the output of this program will be:
```
0.
1.
```
Having thus constructed these objects, we can use any of the normal matrix/vector operators (defined in vector1.h and matrix2.h), for example in the following we create two random matrices with elements between 0. and 1. and do various operations with them:

```
matrix<double> a(3,3);
matrix<double> b(3,3);

for(int i = 0 ; i < 3 ; i ++) {
    for(int j = 0 ; j < 3 ; j ++) {
        a(i,j) = (double)rand()/(double)RAND_MAX;
        b(i,j) = (double)rand()/(double)RAND_MAX;
    }
}

matrix<double> added = a + b; //add the matrices

matrix<double> subtracted = a - b; //subtract b from a

matrix<double> multiplied = a * b; //matrix multiplication

```

where we have introduced the = operator for creating a new matrix and assigning its value to another matrix

we can also do vector/matrix operators, for example, multipliying a random matrix a and a random vector b, and saving the output to a new vector c:

```
matrix<double> a(3,3);
vector1<double> b(3);

for(int i = 0 ; i < 3 ; i ++) {
    for(int j = 0 ; j < 3 ; j ++) {
        a(i,j) = (double)rand()/(double)RAND_MAX;
    }
}

for(int i = 0 ; i < 3 ; i++) {
    b[i] = (double)rand()/(double)RAND_MAX;
}

vector1<double> c = a*b;

```
