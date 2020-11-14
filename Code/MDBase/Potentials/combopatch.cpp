#ifndef COMBOPATCH_CPP
#define COMBOPATCH_CPP

SingPatch::SingPatch(double strr, double disss, double angg) : ComboPatch(1), ang(angg), dis(disss), str(strr)
{
    i1 = new int[2];
    i1[0] = 1;
    i1[1] = 0;
    p = &i1;


    KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(nx, ny, nz, nx, ny, nz, str, dis, ang, 0.75);

    potential_bundle[0] = pot1->clone();


    delete pot1;
    //potential_bundle[0] = pot1;
    //vector1<potentialtheta3D *> pots(1);
}

#endif /* COMBOPATCH_CPP */
