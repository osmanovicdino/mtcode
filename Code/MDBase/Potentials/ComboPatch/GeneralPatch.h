#ifndef GENERALPATCH_H
#define GENERALPATCH_H


struct GeneralPatch : ComboPatch {

    int no_types;
    vector1<int> no_patches_per_type;
    vector1<int> num_per_type;
    vector1<int> total_patches_per_type;
    matrix<double> params;
    matrix<double> orient;
    matrix<int> pot_starters;

    int **i1; //as it is of indeterminate length we have to have a double pointer.

    GeneralPatch(vector1<int>no_patches_per_typee, vector1<int> num_per_typee, matrix<double>&, matrix<double>& );


    ~GeneralPatch() {
        for(int i = 0; i < no_types*(no_types-1)/2 ; i++) {
            delete i1[i];
        }
        delete [] i1;
    }
    
    int num_patches(const int &);
    inline int return_type(int i);
    inline int return_type_patch(int i);
    //inline int mapping_funcion_particles(int i, int j);

    void UpdateIterator(const int &i, const int &j);
    void UpdateIteratorSafe(const int &i, const int &j, int **q);
    int get_total_patches(const int &N);
    void which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj);
    void which_particle(const int &wpi, const int &wpj, int &i, int &j);
    int which_potential(const int &i, const int &j, const int &wpi, const int &wpj);

    void get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);
    
    void CreateFiles();
    GeneralPatch *clone() const
    {
        return new GeneralPatch(*this);
    }
    //the parameter matrix we have must be
};

#include "GeneralPatch.cpp"

#endif /* GENERALPATCH_H */