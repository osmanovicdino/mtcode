// Q lab -> body
// Q^T body -> lab

LangevinNVTR::LangevinNVTR() : LangevinNVT(), im(vector1<double>(1))
{
    //cout << "NVTR constructor called" << endl;
    gammar = 1.0;
    int dimension = this->getgeo().dimension;
    // if(dimension == 1 )
    //     cout << "warning: rotational system cannot be defined with a dimension of 1" << endl;
    matrix<double> rs(1, dimension ); //rotative motion is embedded in D-1 space
    this->setdat(rs);


    matrix<double> res(rs);

    angmom = new matrix<double>;
    *angmom = res;

    orient = new matrix<double>;
    *orient = res;
}

LangevinNVTR::LangevinNVTR(cube &a) : LangevinNVT(a), im(vector1<double>(a.dimension*a.dimension))
{
    //cout << "geo NVTR constructor called" << endl;
    //gammar = 1.0;
    im[0] = 0.25 * 2.0 / 5.0;
    im[4] = 0.25 * 2.0 / 5.0;
    im[8] = 0.25 * 2.0 / 5.0;

    gammar = 1.0;

    Rt = sqrt(2 * gamma * kT / dt);
    Rr = sqrt(2 * gammar * kT / dt);

    int dimension = this->getgeo().dimension;
    if (dimension == 1)
        error("rotational system cannot be defined with a dimension of 1");
    matrix<double> rs(1, dimension); //rotative motion is embedded in D-1 space
    this->setdat(rs);

    matrix<double> res(1,dimension -1);

    angmom = new matrix<double>;
    *angmom = res;

    orient = new matrix<double>;
    *orient = res;

    
}

LangevinNVTR::LangevinNVTR(const LangevinNVTR &a) : LangevinNVT(a), im(vector1<double>(a.dimension * a.dimension))
{
    //cout << "called copy constructor LangevinNVTR" << endl;
    for(int i = 0  ; i < a.im.getsize() ; i++ ) {
        im[i] = (a.im).gpcons(i);   
    }

    gammar = a.gammar;
    Rt = a.Rt;
    Rr = a.Rr;
    dimension = a.dimension;

    this->setdat(a.getdat());

    angmom = new matrix<double>;
    angmom = (*(a.angmom)).clone();

    orient = new matrix<double>;
    orient = (*(a.orient)).clone();
}

LangevinNVTR& LangevinNVTR::operator=(const LangevinNVTR &a) {

    
    delete angmom;
    delete orient;
    LangevinNVT::operator=(a);
    // for (int i = 0; i < a.im.getsize(); i++)
    // {
    //     im[i] = (a.im).gpcons(i);
    // }
    im = a.im;

    gammar = a.gammar;
    Rt = a.Rt;
    Rr = a.Rr;
    dimension = a.dimension;

    this->setdat(a.getdat());

    angmom = new matrix<double>;
    angmom = (*(a.angmom)).clone();
    orient = new matrix<double>;
    orient = (*(a.orient)).clone();

    return *this;
}

LangevinNVTR::~LangevinNVTR() {
    delete angmom;
    delete orient;
}

void LangevinNVTR::setIM(const vector1<double> &I) {
    if(I.getsize() != SQR(dimension)) error("error in forumulation of inertia matrix");
    
    im = I;
}

void LangevinNVTR::setdt(double dtt) {
    dt = dtt;
    d = (gamma * dt / 2.);
    q = (dt) / 2.;
    r = sqrt(0.5 * kT * (gamma) * (m) * (dt));

    c1 = (dt / m);
    c2 = (1.0 / (1.0 + (d)));
    c3 = (1.0 / (1.0 + (d))) * q;
    c4 = (1.0 / (1.0 + (d))) * r;
    c5 = (1 - (d));

    Rt = sqrt(2 * gamma * kT / dt);
    Rr = sqrt(2 * gammar * kT / dt);
}

void LangevinNVTR::measured_temperature(ofstream& myfile) {
    int ds = this->getdimension();
    vector1<double> momav(ds);
    vector1<double> angav(ds);
    for (int i = 0; i < (*mom).getNsafe(); i++)
    {
        for (int i1 = 0; i1 < ds; i1++)
        {
            momav[i1] += ((*mom)(i, i1) * (*mom)(i, i1)) / (2. * m);
            angav[i1] += ((*angmom)(i, i1) * (*angmom)(i, i1)) / (2. * im[0]);
        }
    }
    momav /= ((double)(mom->getNsafe()));
    angav /= ((double)(angmom->getNsafe()));
    
    myfile << momav[0] << ",";
    myfile << momav[1] << ",";
    myfile << momav[2] << "\n";

    myfile << angav[0] << ",";
    myfile << angav[1] << ",";
    myfile << angav[2] << "\n";
}

void LangevinNVTR::measured_temperature(vector1<double> &ret)
{
    int ds = this->getdimension();
    vector1<double> momav(ds);
    vector1<double> angav(ds);
    for (int i = 0; i < (*mom).getNsafe(); i++)
    {
        for (int i1 = 0; i1 < ds; i1++)
        {
            momav[i1] += ((*mom)(i, i1) * (*mom)(i, i1)) / (2. * m);
            angav[i1] += ((*angmom)(i, i1) * (*angmom)(i, i1)) / (2. * im[0]);
        }
    }
    momav /= ((double)(mom->getNsafe()));
    angav /= ((double)(angmom->getNsafe()));

    ret[0] =  momav[0];
    ret[1] = momav[1];
    ret[2] = momav[2];

    ret[3] = angav[0];
    ret[4] = angav[1];
    ret[5] = angav[2];
}

void LangevinNVTR::measured_temperature()
{
    int ds = this->getdimension();
    vector1<double> momav(ds);
    vector1<double> angav(ds);
    for (int i = 0; i < (*mom).getNsafe(); i++)
    {
        for (int i1 = 0; i1 < ds; i1++)
        {
            momav[i1] += ((*mom)(i, i1) * (*mom)(i, i1)) / (2. * m);
            angav[i1] += ((*angmom)(i, i1) * (*angmom)(i, i1)) / (2. * im[0]);
        }
    }
    momav /= ((double)(mom->getNsafe()));
    angav /= ((double)(angmom->getNsafe()));

    cout << momav[0] << ",";
    cout << momav[1] << ",";
    cout << momav[2] << ",";

    cout << angav[0] << ",";
    cout << angav[1] << ",";
    cout << angav[2] << "\n";
}

void LangevinNVTR::initialize(matrix<double> &positions)
{
    this->setdat(positions);

    int NN = this->getN();
    int nc = positions.getncols();


    

    matrix<double> momenta(NN,nc);

    this->setmom(momenta);

    matrix<double> orientations(NN,SQR(nc));
    
    //#pragma omp parallel for
    for(int i = 0 ; i < NN ; i++) {

        double x1 = (double)rand() / (double)(RAND_MAX);
        double x2 = (double)rand() / (double)(RAND_MAX);
        double x3 = (double)rand() / (double)(RAND_MAX);
        double v1 = cos(2 * pi * x2) * sqrt(x3);
        double v2 = sin(2 * pi * x2) * sqrt(x3);
        double v3 = sqrt(1.- x3);

        double r1 = cos(2*pi*x1);
        double r2 = sin(2*pi*x1);

        // // orientations(i, 0) = -(r1 * (1 - SQR(v1))) - r2 * v1 * v2;
        // // orientations(i, 1) = -(r2 * (1 - Power(v1, 2))) + r1 * v1 * v2;
        // // orientations(i, 2) = v1 * v3;
        // // orientations(i, 3) = r1 * v1 * v2 + r2 * (1 - SQR(v2));
        // // orientations(i, 4) = r2 * v1 * v2 - r1 * (1 - SQR(v2));
        // // orientations(i, 5) = v2*v3;
        // // orientations(i, 6) = r1 * v1 * v3 - r2 * v2 * v3;
        // // orientations(i, 7) = r2 * v1 * v3 + r1 * v2 * v3;
        // // orientations(i, 8) = -1 + SQR(v3);

        orientations(i, 0) = -(r1 * (1 - 2 * SQR(v1))) - 2 * r2 * v1 * v2;
        orientations(i, 1) = -(r2 * (1 - 2 * SQR(v1))) + 2 * r1 * v1 * v2;
        orientations(i, 2) = 2 * v1 * v3;
        orientations(i, 3) = 2 * r1 * v1 * v2 + r2 * (1 - 2 * SQR(v2));
        orientations(i, 4) = 2 * r2 * v1 * v2 - r1 * (1 - 2 * SQR(v2));
        orientations(i, 5) = 2 * v2 * v3;
        orientations(i, 6) = 2 * r1 * v1 * v3 - 2 * r2 * v2 * v3;
        orientations(i, 7) = 2 * r2 * v1 * v3 + 2 * r1 * v2 * v3;
        orientations(i, 8) = -1 + 2 * SQR(v3);

        // cout << -orientations(i, 2) * orientations(i, 4) * orientations(i, 6) + orientations(i, 1) * orientations(i, 5) * orientations(i, 6) + orientations(i, 2) * orientations(i, 3) * orientations(i, 7) - orientations(i, 0) * orientations(i, 5) * orientations(i, 7) - orientations(i, 1) * orientations(i, 3) * orientations(i, 8) + orientations(i, 0) * orientations(i, 4) * orientations(i, 8) << endl;
        // cout << orientations.getrowvector(i) << endl;


        // double Psi = 2 * pi * ((double)rand() / (double)RAND_MAX);
        // double Phi = 2 * pi * ((double)rand() / (double)RAND_MAX);
        // double Theta = acos(-1.+2.*((double)rand() / (double)RAND_MAX));
        // orientations(i, 0) = (cos(Psi) * (SQR(cos(Theta)) + SQR(1./tan(Phi))) + SQR(sin(Theta))) * SQR(sin(Phi));
        // orientations(i, 1) = SQR(sin(Theta)) * sin(2 * Phi) * SQR(sin(Psi / 2.)) - cos(Theta) * sin(Psi);
        // orientations(i, 2) = sin(2 * Theta) * sin(Phi) * SQR(sin(Psi / 2.)) + cos(Phi) * sin(Theta) * sin(Psi);
        // orientations(i, 3) = SQR(sin(Theta)) * sin(2 * Phi) * SQR(sin(Psi / 2.)) + cos(Theta) * sin(Psi);
        // orientations(i, 4) = (cos(Psi) + SQR(cos(Theta)) * cos(Psi) * SQR(1./tan(Phi)) + SQR(1./tan(Phi)) * SQR(sin(Theta))) * SQR(sin(Phi));
        // orientations(i, 5) = sin(Theta) * sin(Phi) * (2 * cos(Theta) * (1./tan(Phi)) * SQR(sin(Psi / 2.)) - sin(Psi));
        // orientations(i, 6) = sin(2 * Theta) * sin(Phi) * SQR(sin(Psi / 2.)) - cos(Phi) * sin(Theta) * sin(Psi);
        // orientations(i, 7) = sin(Theta) * sin(Phi) * (2 * cos(Theta) * (1./tan(Phi)) * SQR(sin(Psi / 2.)) + sin(Psi));
        // orientations(i, 8) = SQR(cos(Theta)) + cos(Psi) * SQR(sin(Theta));
        // cout << -orientations(i, 2) * orientations(i, 4) * orientations(i, 6) + orientations(i, 1) * orientations(i, 5) * orientations(i, 6) + orientations(i, 2) * orientations(i, 3) * orientations(i, 7) - orientations(i, 0) * orientations(i, 5) * orientations(i, 7) - orientations(i, 1) * orientations(i, 3) * orientations(i, 8) + orientations(i, 0) * orientations(i, 4) * orientations(i, 8) << endl;


        // cout << orientations.getrowvector(i) << endl;
        // pausel();
    }

    // if (orientations.getNsafe() != this->getN() || orientations.getncols() != SQR(this->getdimension()))
    //     error("incorrect initilization in orientations LangevinNVTR");

    // if (angular_momenta.getNsafe() != this->getN() || angular_momenta.getncols() != this->getdimension())
    //     error("incorrect initilization in orientations LangevinNVTR");

    matrix<double> angular_momenta(NN,3);
    this->setorientation(orientations);
    this->setangularmomenta(angular_momenta);
    // delete orient;
    // orient = new matrix<double>(orientations);

    // delete angmom;
    // angmom = new matrix<double>(angular_momenta);
}

void LangevinNVTR::initialize(matrix<double> &positions, matrix<double> &momenta, matrix<double> &orientations, matrix<double> &angular_momenta) {
    this->setdat(positions);

    
    this->setmom(momenta);

    if (orientations.getNsafe() != this->getN() || orientations.getncols() != SQR(this->getdimension()) )
        error("incorrect initilization in orientations LangevinNVTR");

    if (angular_momenta.getNsafe() != this->getN() || angular_momenta.getncols() != this->getdimension())
        error("incorrect initilization in orientations LangevinNVTR");

    delete orient;
    orient = new matrix<double>(orientations);

    delete angmom;
    angmom = new matrix<double>(angular_momenta);
}

void LangevinNVTR::advancemom_halfstep(matrix<double> &F, matrix<double> &T)  {

    //where F is in LAB FRAME
    //and T is in BODY FRAME
    int np = (*mom).getNsafe();

#pragma omp parallel for schedule(static)
    for (int i = 0; i < np; i++)
    {

            // (mom)->operator()(i, i1) = ((mom)->operator()(i, i1)) + (dt/2.) * F(i, i1) ;
            // (angmom)->operator()(i, i1) = (angmom)->operator()(i, i1) + (dt / 2.) * T(i, i1);
        (*mom)(i, 0) += (dt / 2.) * F(i, 0);
        (*angmom)(i, 0) += (dt / 2.) * T(i, 0);
        (*mom)(i, 1) += (dt / 2.) * F(i, 1);
        (*angmom)(i, 1) += (dt / 2.) * T(i, 1);
        (*mom)(i, 2) += (dt / 2.) * F(i, 2);
        (*angmom)(i, 2) += (dt / 2.) * T(i, 2);
    }


    

    
}

// void LangevinNVTR::advance_pos() {
//     #pragma omp parallel for
//     for (int i = 0; i < (*mom).getNsafe(); i++)
//     {
//         for (int i1 = 0; i1 < dimension; i1++)
//         {
//             (dat)->operator()(i, i1) = ((dat)->operator()(i, i1)) + (dt / m) * mom->operator()(i, i1);
//         }
//     }

//     (*this->geo).correct_position_and_momentum(*dat, *mom);
// }


vector1<double> LangevinNVTR::genfullmat(int i) {
    double ix = im[0];
    double iy = im[4];
    double iz = im[8];

    double jx = angmom->operator()(i, 0);
    double jy = angmom->operator()(i, 1);
    double jz = angmom->operator()(i, 2);

    vector1<double> up(9);

    double tempx = dt * jx / ix;

    double tempy = dt * jy / iy;

    double tempz = dt * jz / iz;

    double fac1 = (1-SQR(tempx/4.))/(1+SQR(tempx/4.));

    double fac2 = tempx/(2.*(1+SQR(tempx/4.)));

    double fac3 = (1 - SQR(tempy / 4.)) / (1 + SQR(tempy / 4.));

    double fac4 = tempy / (2. * (1 + SQR(tempy / 4.)));

    double fac5 = (1 - SQR(tempz / 2.)) / (1 + SQR(tempz / 2.));

    double fac6 = (tempz)/ (1 + SQR(tempz / 2.));

    // up[0] = -(SQR(dt) * SQR(jy)) / (4. * SQR(iy) * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) + (SQR(1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz))));
    // up[1] = -((dt * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) * (1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * jz) / (iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz))))) + (dt * jx * (-(dt * jy * (1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) / (2. * iy * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) - (dt * jy * (1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (2. * iy * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (2. * ix * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix))));
    // up[2] = (SQR(dt) * jx * (1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * jz) / (2. * ix * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) + ((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) * (-(dt * jy * (1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) / (2. * iy * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) - (dt * jy * (1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (2. * iy * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix)));
    // up[3] = -(SQR(dt) * jx * jy * (1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) / (4. * iy * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) + ((1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * ((dt * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jz) / (iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) - (SQR(dt) * jx * jy * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (4. * iy * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)));
    // up[4] = ((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) * ((CUB(dt) * jx * jy * jz) / (4. * iy * SQR(iz) * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) + ((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / ((1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) + (dt * jx * (-(dt * jx * SQR(1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) / (2. * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) - (dt * jy * ((dt * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jz) / (iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) - (SQR(dt) * jx * jy * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (4. * iy * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (2. * iy * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))))) / (2. * ix * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix))));
    // up[5] = -(dt * jx * ((CUB(dt) * jx * jy * jz) / (4. * iy * SQR(iz) * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) + ((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / ((1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (2. * ix * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix)))) + ((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) * (-(dt * jx * SQR(1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) / (2. * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) - (dt * jy * ((dt * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jz) / (iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) - (SQR(dt) * jx * jy * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (4. * iy * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (2. * iy * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))))) / (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix)));
    // up[6] = (dt * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jy * (1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) / (2. * iy * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) + ((1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * ((SQR(dt) * jx * jz) / (2. * SQR(iz) * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) + (dt * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jy * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (2. * iy * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)));
    // up[7] = ((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) * (-(SQR(dt) * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jy * jz) / (2. * iy * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) + (dt * jx * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (2. * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) + (dt * jx * (((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * SQR(1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) / ((1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) - (dt * jy * ((SQR(dt) * jx * jz) / (2. * SQR(iz) * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) + (dt * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jy * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (2. * iy * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (2. * iy * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))))) / (2. * ix * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix))));
    // up[8] = -(dt * jx * (-(SQR(dt) * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jy * jz) / (2. * iy * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) + (dt * jx * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (2. * iz * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (2. * ix * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix)))) + ((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(ix))) * (((1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * SQR(1 - (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) / ((1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * SQR(1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))) - (dt * jy * ((SQR(dt) * jx * jz) / (2. * SQR(iz) * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) + (dt * (1 - (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * jy * (1 - (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))) / (2. * iy * (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(iz))) * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy))) * (1 + (SQR(dt) * SQR(jz)) / (4. * SQR(iz)))))) / (2. * iy * (1 + (SQR(dt) * SQR(jy)) / (16. * SQR(iy)))))) / (1 + (SQR(dt) * SQR(jx)) / (16. * SQR(ix)));
    /*
    double sqrdt = SQR(dt);
    double sqrix = SQR(ix);
    double sqriy = SQR(iy);
    double sqriz = SQR(iz);
    double sqrjx = SQR(jx);
    double sqrjy = SQR(jy);
    double sqrjz = SQR(jz);

    up[0] = -(sqrdt * sqrjy) / (4. * sqriy * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + (SQR(1 - (sqrdt * sqrjy) / (16. * sqriy)) * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (SQR(1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz)));
    up[1] = -((dt * (1 - (sqrdt * sqrjx) / (16. * sqrix)) * (1 - (sqrdt * sqrjy) / (16. * sqriy)) * jz) / (iz * (1 + (sqrdt * sqrjx) / (16. * sqrix)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz)))) + (dt * jx * ((dt * jy * (1 - (sqrdt * sqrjy) / (16. * sqriy))) / (2. * iy * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + (dt * jy * (1 - (sqrdt * sqrjy) / (16. * sqriy)) * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (2. * iy * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (2. * ix * (1 + (sqrdt * sqrjx) / (16. * sqrix)));
    up[2] = (sqrdt * jx * (1 - (sqrdt * sqrjy) / (16. * sqriy)) * jz) / (2. * ix * iz * (1 + (sqrdt * sqrjx) / (16. * sqrix)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) + ((1 - (sqrdt * sqrjx) / (16. * sqrix)) * ((dt * jy * (1 - (sqrdt * sqrjy) / (16. * sqriy))) / (2. * iy * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + (dt * jy * (1 - (sqrdt * sqrjy) / (16. * sqriy)) * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (2. * iy * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (1 + (sqrdt * sqrjx) / (16. * sqrix));
    up[3] = (sqrdt * jx * jy * (1 - (sqrdt * sqrjy) / (16. * sqriy))) / (4. * iy * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + ((1 - (sqrdt * sqrjy) / (16. * sqriy)) * ((dt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jz) / (iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) + (sqrdt * jx * jy * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (4. * iy * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (1 + (sqrdt * sqrjy) / (16. * sqriy));
    up[4] = ((1 - (sqrdt * sqrjx) / (16. * sqrix)) * (-(CUB(dt) * jx * jy * jz) / (4. * iy * sqriz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) + ((1 - (sqrdt * sqrjx) / (16. * sqriz)) * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / ((1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (1 + (sqrdt * sqrjx) / (16. * sqrix)) + (dt * jx * (-(dt * jx * SQR(1 - (sqrdt * sqrjy) / (16. * sqriy))) / (2. * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + (dt * jy * ((dt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jz) / (iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) + (sqrdt * jx * jy * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (4. * iy * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (2. * iy * (1 + (sqrdt * sqrjy) / (16. * sqriy))))) / (2. * ix * (1 + (sqrdt * sqrjx) / (16. * sqrix)));
    up[5] = -(dt * jx * (-(CUB(dt) * jx * jy * jz) / (4. * iy * sqriz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) + ((1 - (sqrdt * sqrjx) / (16. * sqriz)) * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / ((1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (2. * ix * (1 + (sqrdt * sqrjx) / (16. * sqrix))) + ((1 - (sqrdt * sqrjx) / (16. * sqrix)) * (-(dt * jx * SQR(1 - (sqrdt * sqrjy) / (16. * sqriy))) / (2. * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + (dt * jy * ((dt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jz) / (iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) + (sqrdt * jx * jy * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (4. * iy * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (2. * iy * (1 + (sqrdt * sqrjy) / (16. * sqriy))))) / (1 + (sqrdt * sqrjx) / (16. * sqrix));
    up[6] = -(dt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jy * (1 - (sqrdt * sqrjy) / (16. * sqriy))) / (2. * iy * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + ((1 - (sqrdt * sqrjy) / (16. * sqriy)) * ((sqrdt * jx * jz) / (2. * sqriz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) - (dt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jy * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (2. * iy * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (1 + (sqrdt * sqrjy) / (16. * sqriy));
    up[7] = ((1 - (sqrdt * sqrjx) / (16. * sqrix)) * ((sqrdt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jy * jz) / (2. * iy * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) + (dt * jx * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (2. * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (1 + (sqrdt * sqrjx) / (16. * sqrix)) + (dt * jx * (((1 - (sqrdt * sqrjx) / (16. * sqriz)) * SQR(1 - (sqrdt * sqrjy) / (16. * sqriy))) / ((1 + (sqrdt * sqrjx) / (16. * sqriz)) * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + (dt * jy * ((sqrdt * jx * jz) / (2. * sqriz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) - (dt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jy * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (2. * iy * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (2. * iy * (1 + (sqrdt * sqrjy) / (16. * sqriy))))) / (2. * ix * (1 + (sqrdt * sqrjx) / (16. * sqrix)));
    up[8] = -(dt * jx * ((sqrdt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jy * jz) / (2. * iy * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) + (dt * jx * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (2. * iz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (2. * ix * (1 + (sqrdt * sqrjx) / (16. * sqrix))) + ((1 - (sqrdt * sqrjx) / (16. * sqrix)) * (((1 - (sqrdt * sqrjx) / (16. * sqriz)) * SQR(1 - (sqrdt * sqrjy) / (16. * sqriy))) / ((1 + (sqrdt * sqrjx) / (16. * sqriz)) * SQR(1 + (sqrdt * sqrjy) / (16. * sqriy))) + (dt * jy * ((sqrdt * jx * jz) / (2. * sqriz * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))) - (dt * (1 - (sqrdt * sqrjx) / (16. * sqriz)) * jy * (1 - (sqrdt * sqrjz) / (4. * sqriz))) / (2. * iy * (1 + (sqrdt * sqrjx) / (16. * sqriz)) * (1 + (sqrdt * sqrjy) / (16. * sqriy)) * (1 + (sqrdt * sqrjz) / (4. * sqriz))))) / (2. * iy * (1 + (sqrdt * sqrjy) / (16. * sqriy))))) / (1 + (sqrdt * sqrjx) / (16. * sqrix));
    */    

    up[0] = -SQR(fac4) + SQR(fac3) * fac5;
    up[1] = fac3 * (fac2 * fac4 * (1 + fac5) - fac1 * fac6);
    up[2] = fac3 * (fac1 * fac4 * (1 + fac5) + fac2 * fac6);
    up[3] = fac3 * (fac2 * fac4 * (1 + fac5) + fac1 * fac6);
    up[4] = SQR(fac1) * fac5 + SQR(fac2) * (-SQR(fac3) + SQR(fac4) * fac5);
    up[5] = -(fac1 * fac2 * (SQR(fac3) + fac5 - SQR(fac4) * fac5)) + (SQR(fac1) + SQR(fac2)) * fac4 * fac6;
    up[6] = -(fac1 * fac3 * fac4 * (1 + fac5)) + fac2 * fac3 * fac6;
    up[7] = fac1 * fac2 * (SQR(fac3) + fac5 - SQR(fac4) * fac5) + (SQR(fac1) + SQR(fac2)) * fac4 * fac6;
    up[8] = -(SQR(fac2) * fac5) + SQR(fac1) * (SQR(fac3) - SQR(fac4) * fac5);

    return up;
}
// void LangevinNVTR::create_random_forces(matrix<double> &Rt, matrix<double> &Rr) {
//     int Ns = Rt.getNsafe();
//     for(int i = 0 ; i < Ns ; i++) {
//         Rt(i, 0) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
//         Rt(i, 1) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
//         Rt(i, 2) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
//         Rr(i, 0) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
//         Rr(i, 1) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
//         Rr(i, 2) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
//     }
// }

void LangevinNVTR::create_forces_and_torques_sphere(matrix<double> &forcel, matrix<double> &torquel, matrix<double> &Randoms, int starti, int endi)
{ //adds friction and noise
    //int timeseed = int(time(NULL));
    //NO TRANSLATIONAL/ROTATIONAL COUPLING

    //calculate forces and torques in lab frame:

    // int Ns = angmom->getNsafe();
    // int No = orient->getNsafe();
    #pragma omp parallel for schedule(static)
    for (int i = starti; i < endi; i++)
    {
        // srand(timeseed ^ omp_get_thread_num() );
        double qtemp0 = orient->operator()(i, 0);
        double qtemp1 = orient->operator()(i, 1);
        double qtemp2 = orient->operator()(i, 2);
        double qtemp3 = orient->operator()(i, 3);
        double qtemp4 = orient->operator()(i, 4);
        double qtemp5 = orient->operator()(i, 5);
        double qtemp6 = orient->operator()(i, 6);
        double qtemp7 = orient->operator()(i, 7);
        double qtemp8 = orient->operator()(i, 8);

        double plx = mom->operator()(i, 0);
        double ply = mom->operator()(i, 1);
        double plz = mom->operator()(i, 2);

        double llx = angmom->operator()(i, 0);
        double lly = angmom->operator()(i, 1);
        double llz = angmom->operator()(i, 2);

        double pbx = qtemp0 * plx + qtemp1 * ply + qtemp2 * plz;
        double pby = qtemp3 * plx + qtemp4 * ply + qtemp5 * plz;
        double pbz = qtemp6 * plx + qtemp7 * ply + qtemp8 * plz;

        double fbrfx = (-gamma / m) * pbx;
        double fbrfy = (-gamma / m) * pby;
        double fbrfz = (-gamma / m) * pbz;

        double tbrfx = (-gammar / im[0]) * llx;
        double tbrfy = (-gammar / im[4]) * lly;
        double tbrfz = (-gammar / im[8]) * llz;

        double fbrrx = Rt * Randoms(i, 0); // * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
        double fbrry = Rt * Randoms(i, 1); // * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
        double fbrrz = Rt * Randoms(i, 2); // * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
        double tbrrx = Rr * Randoms(i, 3); // * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
        double tbrry = Rr * Randoms(i, 4); // * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
        double tbrrz = Rr * Randoms(i, 5); // * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);


        double fx = forcel(i, 0) + (fbrfx + fbrrx) * qtemp0 + (fbrfy + fbrry) * qtemp3 + (fbrfz + fbrrz) * qtemp6;
        double fy = forcel(i, 1) + (fbrfx + fbrrx) * qtemp1 + (fbrfy + fbrry) * qtemp4 + (fbrfz + fbrrz) * qtemp7;
        double fz = forcel(i, 2) + (fbrfx + fbrrx) * qtemp2 + (fbrfy + fbrry) * qtemp5 + (fbrfz + fbrrz) * qtemp8;


        //double tx = torquel(i, 0) + (tbrfx + tbrrx) * qtemp0 + (tbrfy + tbrry) * qtemp3 + (tbrfz + tbrrz) * qtemp6;
        //double ty = torquel(i, 1) + (tbrfx + tbrrx) * qtemp1 + (tbrfy + tbrry) * qtemp4 + (tbrfz + tbrrz) * qtemp7;
        //double tz = torquel(i, 2) + (tbrfx + tbrrx) * qtemp2 + (tbrfy + tbrry) * qtemp5 + (tbrfz + tbrrz) * qtemp8;
        double tbx = tbrfx + tbrrx + torquel(i, 0) * qtemp0 + torquel(i, 1) * qtemp1 + torquel(i, 2) * qtemp2;
        double tby = tbrfy + tbrry + torquel(i, 0) * qtemp3 + torquel(i, 1) * qtemp4 + torquel(i, 2) * qtemp5;
        double tbz = tbrfz + tbrrz + torquel(i, 0) * qtemp6 + torquel(i, 1) * qtemp7 + torquel(i, 2) * qtemp8;

        // cout << "before" << endl;
        // cout << Rt <<  " " << Rr << endl;
        // cout << forcel(i,0) << " " << torquel(i, 0) << endl;
        // cout << forcel(i, 1) << " " << torquel(i, 1) << endl;
        // cout << forcel(i, 2) << " " << torquel(i, 2) << endl;
        // cout << "after" << endl;
        // cout << fx << " " << tbx << endl;
        // cout << fy << " " << tby << endl;
        // cout << fz << " " << tbz << endl;


        // pausel();

        forcel(i, 0) = fx;
        forcel(i, 1) = fy;
        forcel(i, 2) = fz;

        torquel(i, 0) = tbx;
        torquel(i, 1) = tby;
        torquel(i, 2) = tbz;
        // mom->operator()(i, 0) = mom->operator()(i, 0) + (dt / 2.) * fx;



        // mom->operator()(i, 1) = mom->operator()(i, 1) + (dt / 2.) * fy;
        // mom->operator()(i, 2) = mom->operator()(i, 2) + (dt / 2.) * fz;

        // angmom->operator()(i, 0) = angmom->operator()(i, 0) + (dt / 2.) * tbx;
        // angmom->operator()(i, 1) = angmom->operator()(i, 1) + (dt / 2.) * tby;
        // angmom->operator()(i, 2) = angmom->operator()(i, 2) + (dt / 2.) * tbz;
    }
    //in the following, we assume the centre of resistance is the same as the centre of mass
}

void LangevinNVTR::create_forces_and_torques_sphere(matrix<double> &forcel, matrix<double> &torquel, matrix<double> &Randoms )
{ //adds friction and noise
    //int timeseed = int(time(NULL));
    //NO TRANSLATIONAL/ROTATIONAL COUPLING

    //calculate forces and torques in lab frame:

    int Ns = angmom->getNsafe();
    this->create_forces_and_torques_sphere(forcel, torquel, Randoms, 0,Ns);

    //in the following, we assume the centre of resistance is the same as the centre of mass
}

void LangevinNVTR::rotate() {
    int angn = (angmom)->getNsafe();
    // updates the matrices q and qt;
    #pragma omp parallel for schedule(static)
    for(int i = 0  ; i < angn ; i++) {
        vector1<double> rr = genfullmat(i);

        // cout << rr << endl;
        // pausel();

        double jx = angmom->operator()(i, 0);
        double jy = angmom->operator()(i, 1);
        double jz = angmom->operator()(i, 2);

        double jtemp0 = jx * rr[0] + jy * rr[1] + jz * rr[2];
        double jtemp1 = jx * rr[3] + jy * rr[4] + jz * rr[5];
        double jtemp2 = jx * rr[6] + jy * rr[7] + jz * rr[8];

        // cout << "particle " << i << endl;

        // cout << "Rotate matrix" << endl;

        // cout << rr[0] << " " << rr[1] << " " << rr[2] << endl;
        // cout << rr[3] << " " << rr[4] << " " << rr[5] << endl;
        // cout << rr[6] << " " << rr[7] << " " << rr[8] << endl;

        // cout << endl;

        // cout << "Ang before: " << endl;
        // cout << jx << " " << jy << " " << jz << endl;
        // cout << "Ang after: " << endl;
        // cout << jtemp0 << " " << jtemp1 << " " << jtemp2 << endl;
        // cout << endl;
        

        // double qtemp0 = orient->operator()(i, 0) * rr[0] + orient->operator()(i, 1) * rr[1] + orient->operator()(i, 2) * rr[2];
        // double qtemp1 = orient->operator()(i, 0) * rr[3] + orient->operator()(i, 1) * rr[4] + orient->operator()(i, 2) * rr[5];
        // double qtemp2 = orient->operator()(i, 0) * rr[6] + orient->operator()(i, 1) * rr[7] + orient->operator()(i, 2) * rr[8];
        // double qtemp3 = orient->operator()(i, 3) * rr[0] + orient->operator()(i, 4) * rr[1] + orient->operator()(i, 5) * rr[2];
        // double qtemp4 = orient->operator()(i, 3) * rr[3] + orient->operator()(i, 4) * rr[4] + orient->operator()(i, 5) * rr[5];
        // double qtemp5 = orient->operator()(i, 3) * rr[6] + orient->operator()(i, 4) * rr[7] + orient->operator()(i, 5) * rr[8];
        // double qtemp6 = orient->operator()(i, 6) * rr[0] + orient->operator()(i, 7) * rr[1] + orient->operator()(i, 8) * rr[2];
        // double qtemp7 = orient->operator()(i, 6) * rr[3] + orient->operator()(i, 7) * rr[4] + orient->operator()(i, 8) * rr[5];
        // double qtemp8 = orient->operator()(i, 6) * rr[6] + orient->operator()(i, 7) * rr[7] + orient->operator()(i, 8) * rr[8];

        double qtemp0 = orient->operator()(i,0)* rr[0]+orient->operator()(i,3)* rr[3]+orient->operator()(i,6)* rr[6];

        double qtemp1 = orient->operator()(i, 1) * rr[0] + orient->operator()(i, 4) * rr[3] + orient->operator()(i, 7) * rr[6];

        double qtemp2 = orient->operator()(i, 2) * rr[0] + orient->operator()(i, 5) * rr[3] + orient->operator()(i, 8) * rr[6];

        double qtemp3 = orient->operator()(i, 0) * rr[1] + orient->operator()(i, 3) * rr[4] + orient->operator()(i, 6) * rr[7];

        double qtemp4 = orient->operator()(i, 1) * rr[1] + orient->operator()(i, 4) * rr[4] + orient->operator()(i, 7) * rr[7];

        double qtemp5 = orient->operator()(i, 2) * rr[1] + orient->operator()(i, 5) * rr[4] + orient->operator()(i, 8) * rr[7];

        double qtemp6 = orient->operator()(i, 0) * rr[2] + orient->operator()(i, 3) * rr[5] + orient->operator()(i, 6) * rr[8];

        double qtemp7 = orient->operator()(i, 1) * rr[2] + orient->operator()(i, 4) * rr[5] + orient->operator()(i, 7) * rr[8];

        double qtemp8 = orient->operator()(i, 2) * rr[2] + orient->operator()(i, 5) * rr[5] + orient->operator()(i, 8) * rr[8];

// cout << "Q before: " << endl;
// cout << orient->operator()(i, 0) << " " << orient->operator()(i, 1) << " " << orient->operator()(i, 2) << endl;
// cout << orient->operator()(i, 3) << " " << orient->operator()(i, 4) << " " << orient->operator()(i, 5) << endl;
// cout << orient->operator()(i, 6) << " " << orient->operator()(i, 7) << " " << orient->operator()(i, 8) << endl;
// cout << endl;
// cout << "Q after: " << endl;
// cout << qtemp0 << " " << qtemp1 << " " << qtemp2 << endl;
// cout << qtemp3 << " " << qtemp4 << " " << qtemp5 << endl;
// cout << qtemp6 << " " << qtemp7 << " " << qtemp8 << endl;

// pausel();

        angmom->operator()(i, 0) = jtemp0;
        angmom->operator()(i, 1) = jtemp1;
        angmom->operator()(i, 2) = jtemp2;

        orient->operator()(i, 0) = qtemp0;
        orient->operator()(i, 1) = qtemp1;
        orient->operator()(i, 2) = qtemp2;
        orient->operator()(i, 3) = qtemp3;
        orient->operator()(i, 4) = qtemp4;
        orient->operator()(i, 5) = qtemp5;
        orient->operator()(i, 6) = qtemp6;
        orient->operator()(i, 7) = qtemp7;
        orient->operator()(i, 8) = qtemp8;
    }

}




void LangevinNVTR::calculate_forces_and_torques3D(matrix<int> &pairs, potentialtheta3D &iny, matrix<double> &forces, matrix<double> &torques)
{

int totp =  pairs.getnrows();
#pragma omp parallel for
    for (int i = 0; i < totp; ++i)
    {
        int p1 = pairs(i, 0);
        int p2 = pairs(i, 1);
        //int i1 = pairs(i,2);
        double dis;
        //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
        vector1<double> un(dimension);
        geo.distance_vector(*dat, p1, p2, un, dis);

        //un = i-j
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

        iny.force_and_torque(un, dis, *orient, p1, p2, fx, fy, fz, tix, tiy, tiz, tjx, tjy, tjz);

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
}

    void LangevinNVTR::calculate_forces_and_torques3D(matrix<int> &pairs, vector1<potentialtheta3D*> &iny, matrix<double> &forces, matrix<double> &torques)
    {

        int totp = pairs.getNsafe();
        #pragma omp parallel for
        for (int i = 0; i < totp; ++i)
        {
            int p1 = pairs(i, 0);
            int p2 = pairs(i, 1);
            //int i1 = pairs(i,2);
            double dis;
            //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
            vector1<double> un(dimension);
            geo.distance_vector(*dat, p1, p2, un, dis);

            //un = i-j
            dis = sqrt(dis);

            un /= dis;

            for(int potn = 0 ; potn < iny.getsize() ; potn++) {
                double fx;
                double fy;
                double fz;

                double tix;
                double tiy;
                double tiz;

                double tjx;
                double tjy;
                double tjz;

                iny[potn]->force_and_torque(un, dis, *orient, p1, p2, fx, fy, fz, tix, tiy, tiz, tjx, tjy, tjz);

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
        }
        cout << torques << endl;
        pausel();
    }

    void LangevinNVTR::calculate_forces_and_torques3D(matrix<int> &pairs, ComboPatch &iny, matrix<double> &forces, matrix<double> &torques)
    {
        int totp = pairs.getNsafe();


        #pragma omp parallel for
        for (int i = 0; i < totp; ++i)
        {
            int p1 = pairs(i, 0);
            int p2 = pairs(i, 1);

            if (p2 < p1)
            { //INDICES NEED TO BE SORTED FOR IT TO WORK
                int tp1 = p1;
                p1 = p2;
                p2 = tp1;
            }
            //int i1 = pairs(i,2);
            double dis;
            //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
            vector1<double> un(dimension);
            geo.distance_vector(*dat, p1, p2, un, dis);



            //un = i-j
            dis = sqrt(dis);

            un /= dis;

  


            int **q = new int *;
            
            if (iny.safe)
            {
                iny.UpdateIterator(p1, p2);
                *q = *iny.p;
            }
            else
            {
                iny.UpdateIteratorSafe(p1, p2, q);
            }

            // cout << p1 << " " << p2 << endl;
            // cout << (*q)[0] << endl;
            // cout << (*q)[1] << endl;
            // pausel();
            //iny.UpdateIterator(p1,p2); // this points the patch pointer to the correct particles

            //cout << (*q)[0] << endl;
            for (int tp = 1; tp < (*q)[0] + 1; tp++)
            {
                
                int potn = (*q)[tp];
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

                // cout << p1 << " " << p2 << " " << tp << " " << potn << endl;

                // cout << "forces: " << fx <<" " << fy << " " << fz << endl;
                // cout << dis << endl;
                // cout << un << endl;
                // cout << iny.potential_bundle[potn]->getparameters() << endl;

                // pausel();


                // if((abs(fx)>1E-10 || abs(fy)>1E-10|| abs(fz)> 1E-10) /* &&(p1<2048+50 && p2 >= 2048+50) */ ) {
                // int wp1,wp2;
                // iny.which_patch(p1,p2,potn,wp1,wp2);
                // cout << p1 <<" " << p2 << " " <<potn << " " << wp1 << " " << wp2 << endl;
                // cout << un << endl;
                // cout << dis << endl;
                // cout << fx << " " << fy << " " << fz << endl;
                
                // pausel();
                
                // }
                

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


            delete q;
        }



    }

void LangevinNVTR::calculate_forces_and_torques3D_onlyone(matrix<int> &pairs, vector1<potentialtheta3D *> &iny,  BinaryBindStore &bo, AbstractBindingModel &bm, matrix<double> &forces, matrix<double> &torques)
    {

        //for a given sphere geometry

        int np1 = sqrt(iny.getsize());

        vector1<int> tempbound((this->getN())*np1);

        int depth_of_matrix = 4 ;  //Choose this value to be deep enough such that all values can be stored

        matrix<int> boindices((this->getN())*np1,depth_of_matrix);

        if(SQR(np1) != iny.getsize()) {
            error("number of patches and size of potential bundle incorrect");
        }

        //vector1<int> nump(this->getN(),0); //this is the number of bound particles per particle
        int totp = pairs.getNsafe();
        #pragma omp parallel for
        for (int i = 0; i < totp; ++i)
        {
            int p1 = pairs(i, 0);
            int p2 = pairs(i, 1);

            if(p2<p1) { //INDICES NEED TO BE SORTED FOR IT TO WORK
                int tp1 = p1;
                p1 = p2;
                p2 = tp1;
            }
            //int i1 = pairs(i,2);
            double dis;
            //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
            vector1<double> un(dimension);
            geo.distance_vector(*dat, p1, p2, un, dis);

            //un = i-j
            dis = sqrt(dis);

            un /= dis;

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

            for(int  j = 0 ; j < np1 ; j++) {
                for(int k = 0 ; k < np1 ; k++) {


                    int potn  = np1*j+k;

                    vector1<double> params = iny[potn]->getparameters();

                    int nxb1 = params[0];//iny[potn]->nxb1;
                    int nyb1 = params[1]; //iny[potn]->nyb1;
                    int nzb1 = params[2]; //iny[potn]->nzb1;

                    int nxb2 = params[3]; //iny[potn]->nxb2;
                    int nyb2 = params[4]; //iny[potn]->nyb2;
                    int nzb2 = params[5]; //iny[potn]->nzb2;

                    double disp = params[6];

                    double thetam =  params[8];

                    double nx1 = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
                    double ny1 = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
                    double nz1 = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;

                    double nx2 = nxb2 * gtemp0 + nyb2 * gtemp3 + nzb2 * gtemp6;
                    double ny2 = nxb2 * gtemp1 + nyb2 * gtemp4 + nzb2 * gtemp7;
                    double nz2 = nxb2 * gtemp2 + nyb2 * gtemp5 + nzb2 * gtemp8;

                    double argthetai = -(nx1 * un.gpcons(0) + ny1 * un.gpcons(1) + nz1 * un.gpcons(2));
                    double argthetaj = (nx2 * un.gpcons(0) + ny2 * un.gpcons(1) + nz2 * un.gpcons(2));

                    if (argthetai > cos(thetam) && argthetaj > cos(thetam) && dis < 1.2*disp ) {
                        boindices(p1 * np1 + j, tempbound[p1 * np1 + j]) = p2 * np1 + k;
                        boindices(p2 * np1 + k, tempbound[p2 * np1 + k]) = p1 * np1 + j;

                        tempbound[p1 * np1 + j] += 1;
                        tempbound[p2 * np1 + k] += 1;   
                    }
                }
            }
        }

            //from this get the small world networks

            // vector1<bool> already_accounted((this->getN()) * np1, false);

            // for (int i = 0 ; i < (this->getN()) * np1 ; i++)
            // {
                
            //     if(!already_accounted[i]) {
            //        DFUtil(i,already_accounted,boindices,tempbound);

            //        cout << "\n";
                    
            //     }
            // }

            vector1<int> indexes((this->getN()) * np1);

            vector1<int> nbins = ConnectedComponents(boindices,tempbound,indexes);

            //Now we have the clusters. For each of these clusters, there is so some transition rate from one to another
            for(int i = 0 ; i < nbins.getsize()-1 ; i++) {

                int size_of_cluster = nbins[i+1]-nbins[i];

                if(size_of_cluster ==1) {
                    //do nothing
                    int i1 = indexes[nbins[i]];
                    //isbound[i1]=false;

                    bo.isbound[i1] = false;


                     //not bound to anything.

                }
                else if(size_of_cluster == 2) {
                    //all fine, bindings
                    int ti1 = indexes[nbins[i]];
                    int ti2 = indexes[nbins[i]+1];


                    int i1;
                    int i2;
                    if(ti1 < ti2) {
                        i1 = ti1;
                        i2 = ti2;
                    }
                    else{
                        i1 =  ti2;
                        i2 =  ti1;
                    }


                    bool alreadybound_to_eachother = bo.boundto[i1]==i2 && bo.isbound[i1] && bo.isbound[i2];

                    bool aft;
                    double r = (double)(rand())/(double)(RAND_MAX);
                    bm.doublet(alreadybound_to_eachother, i1, i2, aft, r);

                    if(aft) {
                        bo.boundto[i1] = i2;
                        bo.boundto[i2] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = true;
                    }
                    else{
                        bo.isbound[i1] = false;
                        bo.isbound[i2] = false;
                    }
                    //bool already_bound = prebound(i1,i2);

                }
                else if(size_of_cluster == 3) {
                    //is the cluster fully connected or not?

                    int ti1 = indexes[nbins[i]];
                    int ti2 = indexes[nbins[i]+1];
                    int ti3 = indexes[nbins[i]+2];

                    //SORT THE INDICES (IMPORTANT)

                    int i1;
                    int i2;
                    int i3;

                    sort_triplet(ti1,ti2,ti3,i1,i2,i3);

                    //bool b12,b23,b13;


                    //DETERMINE WHETHER THEY ARE BOUND
                    bool b12 = bo.boundto[i1] == i2 && bo.isbound[i1] && bo.isbound[i2];
                    bool b23 = bo.boundto[i2] == i3 && bo.isbound[i2] && bo.isbound[i3];
                    bool b13 = bo.boundto[i1] == i3 && bo.isbound[i1] && bo.isbound[i3];


                    //DETERMINE THE CONNECTIVENESS OF THE GRAPH
                    //remember, that in order to count as a triplet
                    bool c12 = false;
                    bool c23 = false;
                    bool c13 = false;

                    int nb1 = tempbound[i1];
                    if(nb1 ==1 ) {
                        int tempi = boindices(i1,0);
                        if(tempi == i2) { 
                            c12 = true; 
                            c13 = false;
                            c23 = true; //in order to be a triplet

                        }
                        else if(tempi == i3) { 
                            c13 = true; 
                            c12 = false; 
                            c23 = true;
                            }
                        else error("something weird");

                        //check the other
                        

                    }
                    else if(nb1 ==2) {
                        c12 = true;
                        c13 = true;

                        int nb2 =  tempbound[i2];
                        if(nb1 == 1 ) {
                            c23 = false;
                        }
                        else{
                            c23 = true;
                        }


                    }
                    else{
                        error("error in clustering algorithm");
                    }

                    bool a12;
                    bool a23;
                    bool a13;
                    double r = (double)(rand()) / (double)(RAND_MAX);
                    bm.triplet(b12,b23,b13,c12,c23,c13,i1,i2,i3,a12,a23,a13,r);

                    if(a12) {
                        bo.boundto[i1] = i2;
                        bo.boundto[i2] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = true;
                        bo.isbound[i3] = false;
                    }
                    else if(a23) {
                        bo.boundto[i2] = i3;
                        bo.boundto[i3] = i2;
                        bo.isbound[i1] = false;
                        bo.isbound[i2] = true;
                        bo.isbound[i3] = true;
                    }
                    else if(a13){
                        bo.boundto[i1] = i3;
                        bo.boundto[i3] = i1;
                        bo.isbound[i1] = true;
                        bo.isbound[i2] = false;
                        bo.isbound[i3] = true;
                    }
                    else{
                        bo.isbound[i1] = false;
                        bo.isbound[i2] = false;
                        bo.isbound[i3] = false;
                    }

                }
                else {
                    error("the mythical 4 cluster, get to work, slob!");
                }

            }

            vector1<bool> visited((this->getN()) * np1);

            for (int i = 0; i < (this->getN()) * np1; ++i)
            {
                if(bo.isbound[i] ==true && visited[i]==false ) {
                    visited[bo.boundto[i]] = true;
                    visited[i] = true;
                    int p1 = floor(i/np1); //particle number 1
                    int p2 = floor(bo.boundto[i]/np1); //particle number 2
                    double dis;
                    //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
                    vector1<double> un(dimension);
                    geo.distance_vector(*dat, p1, p2, un, dis);

                    //un = i-j

                    int potn = (i % np1) * np1 + (bo.boundto[i] % np1);
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

                    iny[potn]->force_and_torque(un, dis, *orient, p1, p2, fx, fy, fz, tix, tiy, tiz, tjx, tjy, tjz);

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
                
                else{

                }
            }

            //now we have only the real forces, we no longer need to calculate the forces for the non-bound particles:

            


          
        
    }


    //         matrix<double> forces((*dat).getNsafe(), dimension);
    //         //vec_vec<double> outputs(pairs.getn());
    //         // ofstream myfile;
    //         // myfile.open("forces.csv", ios::out | ios::app);
    //         //cout << pairs.getNsafe() << endl;
    //         //potential * pot = ints(0,1,0).clone();

    // #pragma omp parallel for
    //         for (int i = 0; i < pairs.getNsafe(); ++i)
    //         {
    //             int p1 = pairs.mat[i * 2 + 0];
    //             int p2 = pairs.mat[i * 2 + 1];
    //             //int i1 = pairs(i,2);
    //             double dis;
    //             //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
    //             vector1<double> un(dimension);
    //             geo.distance_vector(*dat, p1, p2, un, dis);

    //             //un = i-j

    //             double f1 = iny.force(sqrt(dis));

    //             for (int j = 0; j < dimension; ++j)
    //             {
    //                 double fac = f1 * un[j] / sqrt(dis);
    //                 (forces).mat[p1 * dimension + j] += fac;
    //                 (forces).mat[p2 * dimension + j] += -fac;
    //             }
    //         }

    //         return forces;
    //     }
    // }





void LangevinNVTR::adv(matrix<int> &pairs) {



}



vector<thetapair> LangevinNVTR::check_arg_thetas_per_pair(matrix<int> &pairs, ComboPatch &iny)
{

    vector<thetapair> edgelist;
    edgelist.reserve(pairs.getnrows());

    unsigned int i;
    int tp = pairs.getNsafe();

    #pragma omp parallel
    {
        vector<thetapair> edgelist_private;
        edgelist_private.reserve(pairs.getnrows());

    #pragma omp for nowait schedule(static)
        for (i = 0; i < tp; ++i)
        {
            int p1 = pairs(i, 0);
            int p2 = pairs(i, 1);
            if (p2 < p1)
            { //INDICES NEED TO BE SORTED FOR IT TO WORK
                int tp1 = p1;
                p1 = p2;
                p2 = tp1;
            }
            //int i1 = pairs(i,2);
            double dis;
            //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
            vector1<double> un(dimension);
            geo.distance_vector(*dat, p1, p2, un, dis);

            //un = i-j

            dis = sqrt(dis);

            if (/*dis < iny.max_check*/ true)
            {
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

                int **q = new int *;
                if (iny.safe)
                {
                    iny.UpdateIterator(p1, p2);
                    *q = *iny.p;
                }
                else
                {
                    iny.UpdateIteratorSafe(p1, p2, q);
                }

                //int **q = iny.p;

                for (int tp = 1; tp < (*q)[0] + 1; tp++)
                {
                    int potn = (*q)[tp];

                    mypot *temppot = iny.potential_bundle[potn];


                    double nxb1 = temppot->nxb1;
                    double nxb2 = temppot->nxb2;
                    double nyb1 = temppot->nyb1;
                    double nyb2 = temppot->nyb2;
                    double nzb1 = temppot->nzb1;
                    double nzb2 = temppot->nzb2;
                    double disp = temppot->interaction_distance;
                    double thetam = temppot->thetam;


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


                        thetapair test;
                        test.thetai = argthetai;
                        test.thetaj = argthetaj;
                        edgelist_private.push_back(test);

                    
                }
                //pausel();
                delete q;
            }
        }
        #pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); i++)
        {
            #pragma omp ordered
            edgelist.insert(edgelist.end(), edgelist_private.begin(), edgelist_private.end());
        }
    }

    return edgelist;
}

vector<patchint> LangevinNVTR::calculate_patch_list(matrix<int> &pairs, ComboPatch &iny)
{
    vector<patchint> edgelist;
    edgelist.reserve(10*pairs.getnrows());
    unsigned int i;
    double max_angle = iny.max_ang;
    double max_dis = iny.max_check + 1.;
    int tp = pairs.getNsafe();

    #pragma omp parallel
    {
        vector<patchint> edgelist_private;
        edgelist_private.reserve(10*pairs.getnrows());

        #pragma omp for nowait schedule(dynamic,100)
        for (i = 0; i < tp; ++i)
        {
            int p1 = pairs(i, 0);
            int p2 = pairs(i, 1);
            if (p2 < p1)
            { //INDICES NEED TO BE SORTED FOR IT TO WORK
                int tp1 = p1;
                p1 = p2;
                p2 = tp1;
            }
            //int i1 = pairs(i,2);
            double dis;
            //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
            vector1<double> un(dimension);
            geo.distance_vector(*dat, p1, p2, un, dis);

            //un = i-j
            

            dis = sqrt(dis);

            if (dis < max_dis)
            {
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

                // vector1<double> qtemp = orient->getrowvector(p1);
                // vector1<double> gtemp = orient->getrowvector(p2);

                // for (int j = 0; j < iny.num_patches(p1) ; j++)
                // {
                //     for (int k = 0; k < iny.num_patches(p2); k++)
                //     {

                //int potn = np1 * j + k;

                int **q = new int *;
                if (iny.safe)
                {
                    iny.UpdateIterator(p1, p2);
                    *q = *iny.p;
                }
                else
                {
                    iny.UpdateIteratorSafe(p1, p2, q);
                }

                //int **q = iny.p;

                for (int tp = 1; tp < (*q)[0] + 1; tp++)
                {
                    int potn = (*q)[tp];

                    mypot *temppot = iny.potential_bundle[potn];

                    double nxb1 = temppot->nxb1;
                    double nxb2 = temppot->nxb2;
                    double nyb1 = temppot->nyb1;
                    double nyb2 = temppot->nyb2;
                    double nzb1 = temppot->nzb1;
                    double nzb2 = temppot->nzb2;
                    double disp = temppot->interaction_distance;
                    double thetam = temppot->thetam;

                    double nx1 = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
                    double ny1 = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
                    double nz1 = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;

                    double nx2 = nxb2 * gtemp0 + nyb2 * gtemp3 + nzb2 * gtemp6;
                    double ny2 = nxb2 * gtemp1 + nyb2 * gtemp4 + nzb2 * gtemp7;
                    double nz2 = nxb2 * gtemp2 + nyb2 * gtemp5 + nzb2 * gtemp8;



                    double argthetai = -(nx1 * dx + ny1 * dy + nz1 * dz);
                    double argthetaj = (nx2 * dx + ny2 * dy + nz2 * dz);

                    if ((argthetai - max_angle) + (argthetaj-max_angle) >-0.4 )
                    {
                        // this is an empirically established formula for how large the total difference betweenn the angles
                        // from the time when they switch on has to be for no new bonds to form after 10 time
                        int wp1, wp2;
                        iny.which_patch(p1, p2, potn, wp1, wp2);

                        patchint test;
                        test.particle_index1 = p1;
                        test.particle_index2 = p2;
                        test.potn = potn;
                        test.patch_index1 = wp1;
                        test.patch_index2 = wp2;
                        edgelist_private.push_back(test);
                    }
                    else{
                        //don't add
                    }
                }
                //pausel();
                delete q;
            }
        }
        #pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); i++)
        {
            #pragma omp ordered
            edgelist.insert(edgelist.end(), edgelist_private.begin(), edgelist_private.end());
        }
    }





    return edgelist;
}

vector<int> adjacency(const vector<patchint> &c1) {
    int n = c1.size();
    vector1<bool> c2(n);
    #pragma omp parallel for 
    for(int i = 1 ; i < n ; i++) {
        bool diff = !(check_same_particle(c1[i],c1[i-1]));
        c2[i] = diff;
    }

    vector<int> indexes;
    indexes.reserve(n/3);
    
    #pragma omp parallel 
    {
    vector<int> indexes_private;
    indexes_private.reserve(n/3);
    #pragma omp for nowait schedule(static)
    for(int i = 0  ; i < n ; i++) {
        if(c2[i]) {
            indexes_private.push_back(i);
        }
    }
    #pragma omp for schedule(static) ordered
    for (int i = 0; i < omp_get_num_threads(); i++)
    {
    #pragma omp ordered
        indexes.insert(indexes.end(), indexes_private.begin(), indexes_private.end());
    }

    }
    #if defined(_OPENMP)
    __gnu_parallel::sort(indexes.begin(), indexes.end());
    #else
        std::sort(indexes.begin(),indexes.end());
    #endif
/* 
    int tn_pairs = c1.size();
    int t_u_pairs = indexes.size();

    vector<patchint> pairs = c1;
    int starti,startj;

    for (int ik = 0; ik < t_u_pairs + 1; ++ik)
    {
        int i, fi; //the start and end indices
        if (ik == 0)
        {
            i = 0;
            fi = indexes[ik];
        }
        else if (ik == t_u_pairs)
        {
            i = indexes[ik - 1];
            fi = tn_pairs;
        }
        else
        {
            i = indexes[ik - 1];
            fi = indexes[ik];
        }

        int p1, p2;
        if(i != startj) {
            cout << i << " " << fi << endl;
            pausel();
        }
        pairs[i].get_particle(p1, p2);
        for (int j = i; j < fi; j++)
            {
            int potn, wp1, wp2;
            pairs[j].get_patch(potn, wp1, wp2);
            cout << p1 << " " << p2 << " " << wp1 << " " << wp2 << " " << potn << endl;
            }
            starti = i;
            startj = fi;
            cout << endl;

        }

        cout << "done" << endl;
        pausel(); */

        return indexes;
    }

void print_single_thetas(vector<thetapair> &a1, string pre, int j)
{
    int n = a1.size();

    stringstream ss;
    ss << j;
    string s1 = "theta";
    string s2 = ss.str();
    string s3 = ".csv";

    ofstream myfile;
    myfile.open((s1 + pre + s2 + s3).c_str());
    for (int i = 0; i < n; i++)
    {
        myfile << a1[i].thetai << "," << a1[i].thetaj << endl;
    }
    myfile.close();
    //pausel();
}

void print_two_different_thetas(vector<thetapair> &a1, vector<thetapair> &a2) {
    int n = a1.size();
    ofstream myfile;
    myfile.open("thetachange.csv");
    for(int i = 0  ; i < n ; i++) {
        myfile << a1[i].thetai - a2[i].thetai << ","  << a1[i].thetaj - a2[i].thetaj << endl;
    
    }
    myfile.close();
    //pausel();
}

void print_two_different_thetas(vector<thetapair> &a1, vector<thetapair> &a2, int j)
{
    int n = a1.size();
    ofstream myfile;
    stringstream ss;
    ss << j;
    string s1 = "thetachange";
    string s2 = ss.str();
    string s3 = ".csv";
    myfile.open((s1+s2+s3).c_str());
    for (int i = 0; i < n; i++)
    {
        myfile << a1[i].thetai - a2[i].thetai << "," << a1[i].thetaj - a2[i].thetaj << endl;
    }
    myfile.close();
    //pausel();
}