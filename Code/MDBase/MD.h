#ifndef MD_H
#define MD_H

#include "../DataStructures/vector1.h"
#include "../DataStructures/vector1.cpp"

#include "../DataStructures/matrix2.h"
#include "../DataStructures/matrix2.cpp"
//#include "matrix.cpp"

#include "../DataStructures/mdpair.h"

// #include "intmatrix.h"
#include "geometry.h"
#include "./Potentials/potential.h"
#include "./Potentials/potential3.h"
#include "./Bindings/GraphAlgorithms.h"

//#include "vec_vec.h"

class MD { //virtual base class of all calculations 3 D
protected:

cube geo; //specify the geometry of the system

matrix<double> *dat; //specifies the positions of all the particles this is the baseline of all molecular dynamics

potential *ints; //interactions between all the elements in the system

int dimension;

public:


MD();
MD(const MD&);
MD& operator=(const MD&);
~MD();

void setdat(const matrix<double>&);

void setinteractions(potential&);
void setgeometry(cube&);

double getcoordinate(int,int);

int getN(); // get the number of particles
int getdimension();
matrix<double> getdat() const;
cube getgeo() const;
potential& getints();
void disvec(int &p1,int&p2,vector1<double>&uv,double&);





double distance(const int&,const int&); //the distance between particles i and j; 
// bool distance_less_than(const int&,const int&, double);
// inline vector1<double> unit_vector(const int&, const int&, double&);
// inline void unit_vector(const int&, const int&, vector1<double>&, double&);
// inline double getcoordinate(int,int);
// void setcoordinate(int,int,double);
//complicated interactions
//void setinteractions();
// double distance(int&, int&); //the distance between particles i and j; 
// bool distance_less_than(const int&,const int&, double);
matrix<int> precalculatepairs(vector<vector<int> > &, matrix<int>&, double);

matrix<int> precalculatepairs(const matrix<int>&, const vector1<int>&, matrix<int> &, double);

vector<mdpair> precalculatepairs_md(const matrix<int> &, const vector1<int> &, matrix<int> &, double);
//matrix<int> GPUprecalculatepairs(vector<vector<int> > &, matrix<int>&, double)

matrix<int>* calculatepairs(matrix<int>&,double); //calculate the pairs over which forces exist

matrix<int> *calculatepairs_parallel(matrix<int> &, double); //calculate the pairs over which forces exist

matrix<int>* calculatepairs(matrix<int>&,vector1<int>&,double); //calculate the pairs over which forces exist, int are the indices of the partial lists 

matrix<int>* calculatepairs(matrix<int>&,vector1<int>&,vector1<int>&,double); //calculate the pairs over which forces exist, int are the indices of the partial lists 

matrix<int>* calculatepairs_sorted(matrix<int>&,double); //calculate the pairs over which forces exist

//virtual matrix<int>* calculatepairs(matrix<int>&);//calculate pairs using the cell list algorithm

//virtual matrix<int>* calculatepairs_global(matrix<int>&,double);//calculate pairs using the cell list algorithm with a global cut off distance

//virtual matrix<int>* calculatepairs_global_onetype(matrix<int>&,double);//calculate pairs using the cell list algorithm with a global cut off distance

//s_matrix<double> calculatepaircorrelation(double,double);

//s_matrix<double> calculateforces_truncated(matrix<int>&,double); //calculate the forces using the pairs as an input, with the forces truncated

//s_matrix<double> calculateforces(matrix<int>&); //calculate the forces using the pairs as an input

matrix<double> calculateforces(matrix<int>&); //calculate the forces using the pairs as an input
matrix<double> calculateforces(matrix<int>&,potential&); //calculate the forces using the pairs as an input
matrix<double> calculateforces_sp(matrix<int>&,potential&,matrix<vector1<double> >&);
matrix<double> calculateforceslist(matrix<int>&,potential&); //calculate the forces using the pairs as an input
matrix<double> calculateforces_threebody(matrix<int>&,potential3&);
matrix<double> calculatestress(matrix<int>&); //calculate the stress;
matrix<double> calculateforces_truncated(matrix<int>&, double); //calculate the forces using the pairs as an input
matrix<double> calculateforces_fast3D(matrix<int>&);

template <class Q>
matrix<double> calculateforces_external(Q&);
//s_matrix<double> calculateforces(matrix<int>&,vector1<int>&,vector1<bool>&,intmatrix&); //calculate the forces using the pairs as an input

//virtual double calculate_energy(vec_vec<int>&);

//virtual double calculate_energy(matrix<int>&);
//virtual void adv(vec_vec<int>&); //advance the system one timestep with the pairs as an input

virtual void adv(matrix<int>&)=0; //advance the system one timestep with the pairs as an input

};



#include "MD.cpp"

#endif