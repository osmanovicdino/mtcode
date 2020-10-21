void set_potential_bundle_bivalent(Condensate &a) {

    double nx = 1;
    double ny = 0;
    double nz = 0;

    KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(nx, ny, nz, nx, ny, nz, 100.0, 2., pi / 4., 0.75);
    KernFrenkelOnePatch2 *pot2 = new KernFrenkelOnePatch2(nx, ny, nz, -nx, ny, nz, 100.0, 2., pi / 4., 0.75);
    KernFrenkelOnePatch2 *pot3 = new KernFrenkelOnePatch2(-nx, ny, nz, nx, ny, nz, 100.0, 2., pi / 4., 0.75);
    KernFrenkelOnePatch2 *pot4 = new KernFrenkelOnePatch2(-nx, ny, nz, -nx, ny, nz, 100.0, 2., pi / 4., 0.75);

    vector1<potentialtheta3D *> pots(4);

    pots[0] = pot1;
    pots[1] = pot2;
    pots[2] = pot3;
    pots[3] = pot4;

    a.set_potential_bundle(pots);

    delete pot1;
    delete pot2;
    delete pot3;
    delete pot4;
}

void set_potential_bundle_tetrahedral(Condensate &a)
{

    //vector1<double> v1(3);

    matrix<double> v(4,3);
    v(0,0) = sqrt(8. / 9.);
    v(0,1) = 0;
    v(0,2) = -1. / 3.;

    //vector1<double> v2(3);

    v(1,0) = -sqrt(2. / 9.);
    v(1,1) = sqrt(2. / 3.);
    v(1,2) = -1. / 3.;
           
    //vector1<double> v3(3);


    v(2,0) = -sqrt(2./9.);
    v(2,1) = -sqrt(2./3.);
    v(2,2) = -1./3.;

    //vector1<double> v4(3);

    v(3,0) = 0;
    v(3,1) = 0;
    v(3,2) = 1.;

    vector1<potentialtheta3D *> pots(16);

    int iter = 0;
    for(int i = 0 ; i < 4 ; i++) {
        for(int j = 0 ; j < 4 ; j++) {
            KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), 100.0, 2., pi / 4., 0.75);

            pots[iter] = pot1;
           // delete pot1;
            iter++;
        }
    }

    a.set_potential_bundle(pots);

    for(int i = 0  ; i < 16 ; i++) {
        delete pots[i];
    }

    // KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(nx, ny, nz, nx, ny, nz, 100.0, 2., pi / 4., 0.75);
    // KernFrenkelOnePatch2 *pot2 = new KernFrenkelOnePatch2(nx, ny, nz, -nx, ny, nz, 100.0, 2., pi / 4., 0.75);
    // KernFrenkelOnePatch2 *pot3 = new KernFrenkelOnePatch2(-nx, ny, nz, nx, ny, nz, 100.0, 2., pi / 4., 0.75);
    // KernFrenkelOnePatch2 *pot4 = new KernFrenkelOnePatch2(-nx, ny, nz, -nx, ny, nz, 100.0, 2., pi / 4., 0.75);

    // vector1<potentialtheta3D *> pots(4);

    // pots[0] = pot1;
    // pots[1] = pot2;
    // pots[2] = pot3;
    // pots[3] = pot4;

    // a.set_potential_bundle(pots);

    // delete pot1;
    // delete pot2;
    // delete pot3;
    // delete pot4;
}