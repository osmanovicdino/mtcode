#ifndef COMBOPATCH_H
#define COMBOPATCH_H

#include "../potentialtheta.h"

typedef KernFrenkelOnePatch2 mypot;

struct ComboPatch {
int **p; //pointer to array where potentials are
bool safe; // is this combo patch parallel safe?
vector1<potentialtheta3D*> potential_bundle; //list of potentials

ComboPatch(int n) : potential_bundle(vector1<potentialtheta3D*>(n)), safe(true) { p = new int*;}

virtual int num_patches(const int&) = 0; //return the number of patches on particle int
virtual void UpdateIterator(const int&, const int&) = 0; //reassign the pointer to point to something else
virtual void UpdateIteratorSafe(const int&, const int&, int**) =0 ;
virtual int get_total_patches(const int &) = 0;
virtual void which_patch(const int&, const int&, const int &, int&, int & ) = 0; //get patch numbers from particle numbers
virtual void which_particle(const int&, const int&, int&, int&)=0; //get particle numbers from patch numbers
virtual int which_potential(const int&, const int&, const int&, const int& )= 0;
virtual void get_params(const int&,const int&, const int&, double&, double&, double&, double&, double&, double& , double&, double&) = 0;
virtual ComboPatch *clone() const = 0;

virtual void CreateFiles();
//~ComboPatch() {delete p;}
};


#include "SingPatch.h"
#include "BivalentPatch.h"
#include "TetrahedralPatch.h"
#include "TetrahedralWithSingle.h"
#include "TetrahedralWithBivalent.h"
#include "TwoTetrahedral.h"
#include "TwoTetrahedralAndSingle.h"
#include "GeneralPatch.h"
#include "combopatch.cpp"
#include "combopatchoutput.cpp"

#endif /* COMBOPATCH_H */
