#ifndef ANALYSISFUNCTIONS_H
#define ANALYSISFUNCTIONS_H

#include "../geometry.h"
#include "../../DataStructures/vector1.h"
#include "../../DataStructures/matrix2.h"

#include "../../Basic/basic.h"
#include "../../DataStructures/mdpair.h"
#include "../Rotational/LangevinR.h"
#include "../Potentials/ComboPatch/combopatch.h"
#include "../Bindings/GraphAlgorithms.h"


vector1<int> calcg1(matrix<double> &data, cube &geo, double max, double dx); //get pair correlation functions

double angd(vector1<double> &angs, int iterator1, int iterator2, double lim); //get the difference in angles (periodicity in 2*pi)


//calculate g1 taking into account the angles
matrix<int> calcg1(matrix<double> &data, vector1<double> &orient, cube &geo, double max, double dx, double dthe);

struct index_test;
vector1<int> distance_graph(matrix<double> &pos, matrix<int> &boxes, cube &geo, double binding_distance, index_test *i1);
//get the growth curve of files in a directory dir
matrix<int> getgrowthcurve(string dir, ComboPatch &iny, int N_Largest);

matrix<int> getgrowthcurve_distance_periodic(string dir, double l, double binding_distance, int N_Largest, index_test *i1);

matrix<double> CreateRandomSample(matrix<double> &pos, int n, double l, int boxes);
//for an initial system with particles at position pos, add n particles at random positions to the sample. l
//is the size of the sample and and boxes is how the initial state should be subdivided for placing particles. 


#include "PairCorrelation.cpp"
#include "AnalyzeGrowth.cpp"
#include "AddParticles.cpp"

#endif /* ANALYSISFUNCTIONS_H */
