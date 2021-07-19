#ifndef GRAPHALGORITHMS_H
#define GRAPHALGORITHMS_H

#include "../../DataStructures/vector1.h"
#include "../../DataStructures/vector1.cpp"

#include "../../DataStructures/matrix2.h"
#include "../../DataStructures/matrix2.cpp"

//iterator for finding connected components and printing them
void DFUtil(int i, vector1<bool> &visited, matrix<int> &adj, vector1<int> &lens);

//iterator for finding connected components and saving them
void DFUtilstore(int i, vector1<bool> &visited, matrix<int> &adj, vector1<int> &lens, vector1<int> &indexes, int &iter, int &len_dem);

/*find connected components of adjacency graph
return the vector of all the demarkers between one component and another
indexes is overwritten with information about the indices in the cluster
connected components will be stored in indexes, where the indices are to be found between the demarkers returned
in the return vector. I.e. for a matrix adj with length lens, the following will print all the indices found together

vector1<int> indexes(200);

vector1<int> nbins = ConnectedComponents(adj,lens,indexes);

for(int i = 0 ; i < nbins.getsize()-1 ; i++) {
    for(int j = nbins[i] ; j < nbins[i+1] ; j++) {
        cout << indexes[j] << " ";
    }
    cout << endl;

}


*/
vector1<int> ConnectedComponents(matrix<int> &adj, vector1<int> &lens, vector1<int> &indexes);





void ConnectedComponentsParallel(matrix<int> &adj, vector1<int> &indexes);

template <typename Y>
void ConnectedComponentsParallel(vector<Y> &adj, vector1<int> &indexes);

template <typename Y>
void ConnectedComponentsParallel(vector<Y> &adj, vector<int> &indexes);

// template <typename Y>
// vector1<int> ConnectedComponentsParallel(vector<Y> &adj, int num_vertices);

template <typename Y>
vector<int> ConnectedComponentsParallel(vector<Y> &adj, int num_vertcies);

bool IndependentEdge(const vector<mdpair> &pairs, const vector1<bool> &bonds);

#include "GraphAlgorithms.cpp"

#endif