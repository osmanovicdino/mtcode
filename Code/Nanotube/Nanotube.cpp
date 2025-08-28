#ifndef NANOTUBE_CPP
#define NANOTUBE_CPP


struct mypoin {

vector1<double> operator()(double r, double theta, double phi) {
    //{r Cos[\[Phi]] Sin[\[Theta]], r Sin[\[Theta]] Sin[\[Phi]], r Cos[\[Theta]]}
    vector1<double> a(3);
    a[0] = r*cos(phi)*sin(theta);
    a[1] = r*sin(phi)*sin(theta);
    a[2] = r*cos(theta);
    return a;


}

};

GeneralPatch CreateStoppers(double int1 ,  double size, double range, double ang, int m1, int m2, int n) {

vector1<int> vec1(3);
vec1[0] = 4;
vec1[1] = 2;
vec1[2] = 2;

vector1<int> numb(3);

numb[0] = m1;
numb[1] = m2;
numb[2] = n;

int tot = 4*4 + 4*2 + 4*2 + 2*2 + 2*2 + 2*2;
matrix<double> params(tot, 3);



int iter = 0;
for (int i = 0; i < 4; i++) 
{
    for (int j = 0; j < 4; j++)
    {
        params(iter, 0) = 0.; //junction-junction interaction
        params(iter, 1) = range * size;
        params(iter, 2) = ang;
        iter++;
    }
}

for (int i = 0; i < 4; i++) 
{
    for (int j = 0; j < 2; j++)
    {

        params(iter, 0) = int1; //junction-tube interaction
        params(iter, 1) = range * size;
        params(iter, 2) = ang;

        iter++;
    }
}

for (int i = 0; i < 4; i++) 
{
    for (int j = 0; j < 2; j++)
    {

            params(iter, 0) = int1; //junction-stopper interaction
            params(iter, 1) = range * size;
            params(iter, 2) = ang; // slightly smaller aperture
        
        iter++;
    }
}
for (int i = 0; i < 2; i++) // tube-tube interaction
{
    for (int j = 0; j < 2; j++)
    {

            params(iter, 0) = int1;
            params(iter, 1) = range * size;
            params(iter, 2) = ang;
        

        iter++;
    }
}

for (int i = 0; i < 2; i++) //tube-stopper interaction
{
    for (int j = 0; j < 2; j++)
    {


        params(iter, 0) = int1;
        params(iter, 1) = range * size; // destroy the gas phase
        params(iter, 2) = ang;
        iter++;
    }

}

for (int i = 0; i < 2; i++) // stopper-stopper interaction
{
    for (int j = 0; j < 2; j++)
    {

        params(iter, 0) = 0.0;
        params(iter, 1) = range * size;
        params(iter, 2) = ang;

        iter++;
    }
}

matrix<double> orient(4+2+2, 3);


orient(0, 0) = -1.0;
orient(0, 1) = 0.;
orient(0, 2) = 0.;

orient(1, 0) = 1.0;
orient(1, 1) = 0.;
orient(1, 2) = 0.;

orient(2, 0) = 0.;
orient(2, 1) = 1.;
orient(2, 2) = 0.;

orient(3, 0) = 0.;
orient(3, 1) = -1.;
orient(3, 2) = 0.;

orient(4, 0) = 1.0;
orient(4, 1) = 0.;
orient(4, 2) = 0.;

orient(5, 0) = -1.;
orient(5, 1) = 0.;
orient(5, 2) = 0.;

orient(6, 0) = 1.;
orient(6, 1) = 0.;
orient(6, 2) = 0.;

orient(7, 0) = 0.;
orient(7, 1) = 1.;
orient(7, 2) = 0.;

// orient(8, 0) = nx5;
// orient(8, 1) = ny5;
// orient(8, 2) = nz5;

// orient(9, 0) = nx6;
// orient(9, 1) = ny6;
// orient(9, 2) = nz6;

// orient(10, 0) = nx7;
// orient(10, 1) = ny7;
// orient(10, 2) = nz7;

// orient(11, 0) = nx8;
// orient(11, 1) = ny8;
// orient(11, 2) = nz8;


// orient(11, 0) = 0.;
// orient(11, 1) = 0.;
// orient(11, 2) = 1.;

GeneralPatch c(vec1, numb, params, orient);
return c;
}

GeneralPatch CreateHexatic(int N, double str1, double rang1, double ang1,
                                double str2, double rang2, double ang2,
                                double dphi, double dtheta,
                                double baseangle
                                )
{

    int nop = 12;
    //pots = new ComboPatch;
    vector1<int> vec1(1);
    vec1[0] = nop;

    vector1<int> numb(1);

    numb[0] = N;

    int tot = nop * nop;
    matrix<double> params(tot, 3);

    int iter = 0;
    for (int i = 0; i < nop; i++)
    {
        for (int j = 0; j < nop; j++)
        {
            bool cond1 = i < 6 && j < 6;
            bool cond2 = i >= 6 && j >= 6 && i < 12 && j < 12;
            bool cond3 = i >= 12 && j >= 12;

            if (cond1 || cond2)
            {
                params(iter, 0) = str1;
                params(iter, 1) = rang1;
                params(iter, 2) = ang1;
                iter++;
            }
            else if (cond3)
            {
                params(iter, 0) = 30.0;
                params(iter, 1) = 1.2;
                params(iter, 2) = 0.4;
                iter++;
            }
            else
            {
                params(iter, 0) = str2;
                params(iter, 1) = rang2;
                params(iter, 2) = ang2;
                iter++;
            }
        }
    }

    outfunc(params, "myparameters");

    matrix<double> orient(nop, 3);

    mypoin B;
    // orient(0, 0) = 0.;
    // orient(0, 1) = 0.;
    // orient(0, 2) = 1.;

    // orient(1, 0) = 1;
    // orient(1, 1) = 0.;
    // orient(1, 2) = 0.;

    // orient(2, 0) = -0.5;
    // orient(2, 1) = sqrt(3)/2;
    // orient(2, 2) = 0.;

    // orient(3, 0) = 0.;
    // orient(3, 1) = 0.;
    // orient(3, 2) = -1.;

    // double dphi = 0.2;
    // double dtheta = 0.1;
    vector1<double> p1 = B(1, dtheta + pid / 2., -dphi);
    vector1<double> p2 = B(1, pid / 2., -dphi);
    vector1<double> p3 = B(1, -dtheta + pid / 2., -dphi);

    vector1<double> p4 = B(1, dtheta + pid / 2., baseangle + dphi);
    vector1<double> p5 = B(1, pid / 2., baseangle + dphi);
    vector1<double> p6 = B(1, -dtheta + pid / 2., baseangle + dphi);

    vector1<double> p7 = B(1, dtheta + pid / 2., dphi);
    vector1<double> p8 = B(1, pid / 2., dphi);
    vector1<double> p9 = B(1, -dtheta + pid / 2., dphi);

    vector1<double> p10 = B(1, dtheta + pid / 2., baseangle - dphi);
    vector1<double> p11 = B(1, pid / 2., baseangle - dphi);
    vector1<double> p12 = B(1, -dtheta + pid / 2., baseangle - dphi);

    // vector1<double> p13(3);
    // vector1<double> p14(3);

    // double p1[3] = {((1 + sqrt(3)) * cos(pid / 18.)) / (2. * sqrt(2)), -((1 + sqrt(3)) * sin(pid / 18.)) / (2. * sqrt(2)), -(-1 + sqrt(3)) / (2. * sqrt(2))};

    // double p2[3] = {cos(pid / 18.), -sin(pid / 18.), 0};

    // double p3[3] = {((1 + sqrt(3)) * cos(pid / 18.)) / (2. * sqrt(2)), -((1 + sqrt(3)) * sin(pid / 18.)) / (2. * sqrt(2)), (-1 + sqrt(3)) / (2. * sqrt(2))};

    // double p4[3] = {-((1 + sqrt(3)) * sin((2 * pid) / 9.)) / (2. * sqrt(2)), ((1 + sqrt(3)) * cos((2 * pid) / 9.)) / (2. * sqrt(2)), -(-1 + sqrt(3)) / (2. * sqrt(2))};

    // double p5[3] = {-sin((2 * pid) / 9.), cos((2 * pid) / 9.), 0};

    // double p6[3] = {-((1 + sqrt(3)) * sin((2 * pid) / 9.)) / (2. * sqrt(2)), ((1 + sqrt(3)) * cos((2 * pid) / 9.)) / (2. * sqrt(2)), (-1 + sqrt(3)) / (2. * sqrt(2))};

    // double p7[3] = {((1 + sqrt(3)) * cos(pid / 18.)) / (2. * sqrt(2)), ((1 + sqrt(3)) * sin(pid / 18.)) / (2. * sqrt(2)), -(-1 + sqrt(3)) / (2. * sqrt(2))};

    // double p8[3] = {cos(pid / 18.), sin(pid / 18.), 0};

    // double p9[3] = {((1 + sqrt(3)) * cos(pid / 18.)) / (2. * sqrt(2)), ((1 + sqrt(3)) * sin(pid / 18.)) / (2. * sqrt(2)), (-1 + sqrt(3)) / (2. * sqrt(2))};

    // double p10[3] = {-((1 + sqrt(3)) * sin(pid / 9.)) / (2. * sqrt(2)), ((1 + sqrt(3)) * cos(pid / 9.)) / (2. * sqrt(2)), -(-1 + sqrt(3)) / (2. * sqrt(2))};

    // double p11[3] = {-sin(pid / 9.), cos(pid / 9.), 0};

    // double p12[3] = {-((1 + sqrt(3)) * sin(pid / 9.)) / (2. * sqrt(2)), ((1 + sqrt(3)) * cos(pid / 9.)) / (2. * sqrt(2)), (-1 + sqrt(3)) / (2. * sqrt(2))};

    orient(0, 0) = p1[0];
    orient(0, 1) = p1[1];
    orient(0, 2) = p1[2];

    orient(1, 0) = p2[0];
    orient(1, 1) = p2[1];
    orient(1, 2) = p2[2];

    orient(2, 0) = p3[0];
    orient(2, 1) = p3[1];
    orient(2, 2) = p3[2];

    orient(3, 0) = p4[0];
    orient(3, 1) = p4[1];
    orient(3, 2) = p4[2];

    orient(4, 0) = p5[0];
    orient(4, 1) = p5[1];
    orient(4, 2) = p5[2];

    orient(5, 0) = p6[0];
    orient(5, 1) = p6[1];
    orient(5, 2) = p6[2];

    orient(6, 0) = p7[0];
    orient(6, 1) = p7[1];
    orient(6, 2) = p7[2];

    orient(7, 0) = p8[0];
    orient(7, 1) = p8[1];
    orient(7, 2) = p8[2];

    orient(8, 0) = p9[0];
    orient(8, 1) = p9[1];
    orient(8, 2) = p9[2];

    orient(9, 0) = p10[0];
    orient(9, 1) = p10[1];
    orient(9, 2) = p10[2];

    orient(10, 0) = p11[0];
    orient(10, 1) = p11[1];
    orient(10, 2) = p11[2];

    orient(11, 0) = p12[0];
    orient(11, 1) = p12[1];
    orient(11, 2) = p12[2];

    // orient(12,0) = 0.;
    // orient(12,1) = 0.;
    // orient(12,2) = 1.0;

    // orient(13, 0) = 0.;
    // orient(13, 1) = 0.;
    // orient(13, 2) = -1.0;

    GeneralPatch c(vec1, numb, params, orient, true);
    return c;
}

NanotubeAssembly::NanotubeAssembly(double rmax, int N) {


    ll =  2*rmax;

    num = floor(ll/ 4.);

    // cout << "beginning to set" << endl;

    obj = new LangevinNVTR;


    double deltaG =  20.;
    double angle = 0.5;

    int tot = 4 * 4 + 4 * 2 + 2 * 2;
    matrix<double> params(tot, 3);

    for (int i = 0; i < tot; i++)
    {
        params(i, 0) = deltaG;
        params(i, 1) = 1.4;
        params(i, 2) = angle;
    }

    TetrahedralWithBivalent c2(params, N, N);


    
    this->setpots(c2);



    spherical_confinement_3D conf2(rmax,1.0,ll/2.);
    // conf2.rmax = rmax/2.;
    // conf.v = 1.0;
    // conf.shft  = ll/2.;
    conf = conf2;
    // conf.rmax = rmax/2.;
    // conf.v = 1.0;
    // conf.shft =  ll/2.;
    // conf = conf2;

    vector1<bool> pb(3,  false);
    cube geo(ll, pb, 3);

    matrix<double> dat(N, 3);

    int pp = floor(ll - 1.);


    vector<double> possible_pos_x;
    vector<double> possible_pos_y;
    vector<double> possible_pos_z;


    myrmax = 0.4*rmax;


    for (int i = 0; i < pp; i++)
    {
        for (int j = 0; j < pp; j++)
        {
            for (int k = 0; k < pp; k++)
            {
                double x = 0.5 + i;
                double y = 0.5 + j;
                double z = 0.5 + k;
                if (SQR(x - ll / 2.) + SQR(y - ll / 2.) + SQR(z - ll / 2.) < SQR(0.9 * myrmax) )
                {
                    //cout << i <<  " " << j << " " << k << endl;
                    possible_pos_x.push_back(x);
                    possible_pos_y.push_back(y);
                    possible_pos_z.push_back(z);
                }
            }
        }
    }




    for (int i = 0; i < N; i++)
    {
        int randint = rand() % (possible_pos_x.size());
        dat(i, 0) = possible_pos_x[randint];
        dat(i, 1) = possible_pos_y[randint];
        dat(i, 2) = possible_pos_z[randint];

        possible_pos_x.erase(possible_pos_x.begin() + randint);
        possible_pos_y.erase(possible_pos_y.begin() + randint);
        possible_pos_z.erase(possible_pos_z.begin() + randint);
    }


    // cout << "beginning" << endl;

    LangevinNVTR b(geo);

    // cout << "created" << endl;
    outfunc(dat,"init");

    b.initialize(dat);

    double kT = 1.0;
    double dt = 0.005;
    b.setdt(dt);

    double viscosity = 1.0;
    double hdradius = 0.5;

    b.setgamma(6. * pi * viscosity * hdradius);
    b.setgammar(8. * pi * viscosity * hdradius * hdradius * hdradius);

    b.setkT(kT);

    *obj = b;

}

NanotubeAssembly::NanotubeAssembly(double llx, int N, bool cc)
{

    ll = llx;

    num = floor(ll / 4.);

    // cout << "beginning to set" << endl;

    obj = new LangevinNVTR;

    double deltaG = 20.;
    double angle = 0.5;

    int tot = 4 * 4 + 4 * 2 + 2 * 2;
    matrix<double> params(tot, 3);

    for (int i = 0; i < tot; i++)
    {
        params(i, 0) = deltaG;
        params(i, 1) = 1.4;
        params(i, 2) = angle;
    }

    TetrahedralWithBivalent c2(params, N, N);

    this->setpots(c2);

    spherical_confinement_3D conf2(ll, 1.0, ll / 2.);
    // conf2.rmax = rmax/2.;
    // conf.v = 1.0;
    // conf.shft  = ll/2.;
    conf = conf2;
    // conf.rmax = rmax/2.;
    // conf.v = 1.0;
    // conf.shft =  ll/2.;
    // conf = conf2;

    vector1<bool> pb(3, false);
    cube geo(ll, pb, 3);

    matrix<double> dat(N, 3);

    int pp = floor(ll - 1.);

    vector<double> possible_pos_x;
    vector<double> possible_pos_y;
    vector<double> possible_pos_z;


    for (int i = 0; i < N; i++)
    { //we do not care about overlaps
        
        dat(i, 0) = ll * (double)rand() / (double)RAND_MAX;
        dat(i, 1) = ll * (double)rand() / (double)RAND_MAX;
        dat(i, 2) = ll * (double)rand() / (double)RAND_MAX;
    }

    // cout << "beginning" << endl;

    LangevinNVTR b(geo);

    // cout << "created" << endl;
   // outfunc(dat, "init");

    b.initialize(dat);

    double kT = 1.0;
    double dt = 0.005;
    b.setdt(dt);

    double viscosity = 1.0;
    double hdradius = 0.5;

    b.setgamma(6. * pi * viscosity * hdradius);
    b.setgammar(8. * pi * viscosity * hdradius * hdradius * hdradius);

    b.setkT(kT);

    *obj = b;
}

void NanotubeAssembly::add_particle2() {

    double r1 = ((double)rand() / (double)(RAND_MAX));
    double r2 = ((double)rand() / (double)(RAND_MAX));
    double r3 = ((double)rand() / (double)(RAND_MAX));
    double x1 = ll / 2. + (0.9 * myrmax / 2.) * (2.*r1-1.);
    double y1 = ll / 2. + (0.9 * myrmax / 2.) * (2.*r2-1.);
    double z1 = ll / 2. + (0.9 * myrmax / 2.) * (2.*r3-1.);

    vector1<double> myvec1(3);
    myvec1[0] = x1;
    myvec1[1] = y1;
    myvec1[2] = z1;

    matrix<double> olddat = obj->getdat();
    matrix<double> oldmom = obj->getmom();
    matrix<double> oldorient = obj->getorientation();
    matrix<double> oldamom = obj->getangmom();

    vector1<double> mydistances(obj->getN());
    for (int i = 0; i < obj->getN(); i++) {
        double x2 = olddat(i, 0);
        double y2 = olddat(i, 1);
        double z2 = olddat(i, 2);
        vector1<double> myvec2(3);
        myvec2[0] = x2;
        myvec2[1] = y2;
        myvec2[2] = z2;
        mydistances[i] = obj->getgeo().distance(myvec1, myvec2);
    }

    double mymin = minval(mydistances);



    int iter = 0;
    while(mymin<1.0) {
        r1 = ((double)rand() / (double)(RAND_MAX));
        r2 = ((double)rand() / (double)(RAND_MAX));
        r3 = ((double)rand() / (double)(RAND_MAX));
        x1 = ll / 2. + (0.9 * myrmax / 2.) * (2*r1-1.);
        y1 = ll / 2. + (0.9 * myrmax / 2.) * (2*r2-1.);
        z1 = ll / 2. + (0.9 * myrmax / 2.) * (2*r3-1.);

        myvec1[0] = x1;
        myvec1[1] = y1;
        myvec1[2] = z1;

        for (int i = 0; i < obj->getN(); i++)
        {
            double x2 = olddat(i, 0);
            double y2 = olddat(i, 1);
            double z2 = olddat(i, 2);
            vector1<double> myvec2(3);
            myvec2[0] = x2;
            myvec2[1] = y2;
            myvec2[2] = z2;
            mydistances[i] = obj->getgeo().distance(myvec1, myvec2);
        }
        mymin = minval(mydistances);
    }
    int NN = obj->getN();
    matrix<double> newdat(NN+1,3);
    for(int i = 0 ; i < NN ; i++) {
        newdat(i, 0) = olddat(i,0);
        newdat(i, 1) = olddat(i, 1);
        newdat(i, 2) = olddat(i, 2);
    }

    newdat(NN, 0) = x1;
    newdat(NN, 1) = y1;
    newdat(NN, 2) = z1;
    oldmom.resize_keep(NN+1,3);
    oldorient.resize_keep(NN + 1, 9);
    oldamom.resize_keep(NN + 1, 3);

    oldorient(NN, 0) = 1.0;
    oldorient(NN, 4) = 1.0;
    oldorient(NN, 8) = 1.0;
    obj->setdat(newdat);
    obj->setmom(oldmom);
    obj->setorientation(oldorient);
    obj->setangularmomenta(oldamom);
    // if (SQR(x - ll / 2.) + SQR(y - ll / 2.) + SQR(z - ll / 2.) < SQR(0.9 * myrmax / 2.))
}


/*
void NanotubeAssembly::add_particle42(int which)
{

    double r1 = ((double)rand() / (double)(RAND_MAX));
    double r2 = ((double)rand() / (double)(RAND_MAX));
    double r3 = ((double)rand() / (double)(RAND_MAX));
    double x1 = ll / 2. + (0.9 * myrmax / 2.) * (2. * r1 - 1.);
    double y1 = ll / 2. + (0.9 * myrmax / 2.) * (2. * r2 - 1.);
    double z1 = ll / 2. + (0.9 * myrmax / 2.) * (2. * r3 - 1.);

    if(which==0) {
        //add a 4

        vector1<double> myvec1(3);
        myvec1[0] = x1;
        myvec1[1] = y1;
        myvec1[2] = z1;

        matrix<double> olddat = obj->getdat();
        matrix<double> oldmom = obj->getmom();
        matrix<double> oldorient = obj->getorientation();
        matrix<double> oldamom = obj->getangmom();

        vector1<double> mydistances(obj->getN());
        for (int i = 0; i < obj->getN(); i++)
        {
            double x2 = olddat(i, 0);
            double y2 = olddat(i, 1);
            double z2 = olddat(i, 2);
            vector1<double> myvec2(3);
            myvec2[0] = x2;
            myvec2[1] = y2;
            myvec2[2] = z2;
            mydistances[i] = obj->getgeo().distance(myvec1, myvec2);
        }

        double mymin = minval(mydistances);

        int iter = 0;
        while (mymin < 1.0)
        {
            r1 = ((double)rand() / (double)(RAND_MAX));
            r2 = ((double)rand() / (double)(RAND_MAX));
            r3 = ((double)rand() / (double)(RAND_MAX));
            x1 = ll / 2. + (0.9 * myrmax / 2.) * (2 * r1 - 1.);
            y1 = ll / 2. + (0.9 * myrmax / 2.) * (2 * r2 - 1.);
            z1 = ll / 2. + (0.9 * myrmax / 2.) * (2 * r3 - 1.);

            myvec1[0] = x1;
            myvec1[1] = y1;
            myvec1[2] = z1;

            for (int i = 0; i < obj->getN(); i++)
            {
                double x2 = olddat(i, 0);
                double y2 = olddat(i, 1);
                double z2 = olddat(i, 2);
                vector1<double> myvec2(3);
                myvec2[0] = x2;
                myvec2[1] = y2;
                myvec2[2] = z2;
                mydistances[i] = obj->getgeo().distance(myvec1, myvec2);
            }
            mymin = minval(mydistances);
        }
        int NN = obj->getN();
        matrix<double> newdat(NN + 1, 3);
        for (int i = 0; i < NN; i++)
        {
            newdat(i, 0) = olddat(i, 0);
            newdat(i, 1) = olddat(i, 1);
            newdat(i, 2) = olddat(i, 2);
        }

        newdat(NN, 0) = x1;
        newdat(NN, 1) = y1;
        newdat(NN, 2) = z1;
        oldmom.resize_keep(NN + 1, 3);
        oldorient.resize_keep(NN + 1, 9);
        oldamom.resize_keep(NN + 1, 3);

        oldorient(NN, 0) = 1.0;
        oldorient(NN, 4) = 1.0;
        oldorient(NN, 8) = 1.0;

        int Nt = pots->nt;

        vector1<double> v1 = newdat.getrowvector(Nt);
        vector1<double> v2 = oldmom.getrowvector(Nt);
        vector1<double> v3 = oldorient.getrowvector(Nt);
        vector1<double> v4 = oldamom.getrowvector(Nt);

        vector1<double> u1 = newdat.getrowvector(NN);
        vector1<double> u2 = oldmom.getrowvector(NN);
        vector1<double> u3 = oldorient.getrowvector(NN);
        vector1<double> u4 = oldamom.getrowvector(NN);

        //swap the rows

        // cout << NN + 1 << endl;
        // cout << Nt << endl;

        // for(int i= 2048 ; i < NN+1 ; i++)
        // cout << newdat.getrowvector(i) << endl;
        // cout << endl;

        newdat.setrow(u1, Nt);
        oldmom.setrow(u2, Nt);
        oldorient.setrow(u3, Nt);
        oldamom.setrow(u4, Nt);

        newdat.setrow(v1, NN);
        oldmom.setrow(v2, NN);
        oldorient.setrow(v3, NN);
        oldamom.setrow(v4, NN);

        (*pots).nt = Nt+1; //set the new pot

        // for (int i = 2048; i < NN + 1; i++)
        //     cout << newdat.getrowvector(i) << endl;
        // cout << endl;

        // pausel();

        obj->setdat(newdat);
        obj->setmom(oldmom);
        obj->setorientation(oldorient);
        obj->setangularmomenta(oldamom);
    }
    else if(which ==1) {
        //add a 2

        vector1<double> myvec1(3);
        myvec1[0] = x1;
        myvec1[1] = y1;
        myvec1[2] = z1;

        matrix<double> olddat = obj->getdat();
        matrix<double> oldmom = obj->getmom();
        matrix<double> oldorient = obj->getorientation();
        matrix<double> oldamom = obj->getangmom();

        vector1<double> mydistances(obj->getN());
        for (int i = 0; i < obj->getN(); i++)
        {
            double x2 = olddat(i, 0);
            double y2 = olddat(i, 1);
            double z2 = olddat(i, 2);
            vector1<double> myvec2(3);
            myvec2[0] = x2;
            myvec2[1] = y2;
            myvec2[2] = z2;
            mydistances[i] = obj->getgeo().distance(myvec1, myvec2);
        }

        double mymin = minval(mydistances);

        int iter = 0;
        while (mymin < 1.0)
        {
            r1 = ((double)rand() / (double)(RAND_MAX));
            r2 = ((double)rand() / (double)(RAND_MAX));
            r3 = ((double)rand() / (double)(RAND_MAX));
            x1 = ll / 2. + (0.9 * myrmax / 2.) * (2 * r1 - 1.);
            y1 = ll / 2. + (0.9 * myrmax / 2.) * (2 * r2 - 1.);
            z1 = ll / 2. + (0.9 * myrmax / 2.) * (2 * r3 - 1.);

            myvec1[0] = x1;
            myvec1[1] = y1;
            myvec1[2] = z1;

            for (int i = 0; i < obj->getN(); i++)
            {
                double x2 = olddat(i, 0);
                double y2 = olddat(i, 1);
                double z2 = olddat(i, 2);
                vector1<double> myvec2(3);
                myvec2[0] = x2;
                myvec2[1] = y2;
                myvec2[2] = z2;
                mydistances[i] = obj->getgeo().distance(myvec1, myvec2);
            }
            mymin = minval(mydistances);
        }
        int NN = obj->getN();
        matrix<double> newdat(NN + 1, 3);
        for (int i = 0; i < NN; i++)
        {
            newdat(i, 0) = olddat(i, 0);
            newdat(i, 1) = olddat(i, 1);
            newdat(i, 2) = olddat(i, 2);
        }

        newdat(NN, 0) = x1;
        newdat(NN, 1) = y1;
        newdat(NN, 2) = z1;
        oldmom.resize_keep(NN + 1, 3);
        oldorient.resize_keep(NN + 1, 9);
        oldamom.resize_keep(NN + 1, 3);

        oldorient(NN, 0) = 1.0;
        oldorient(NN, 4) = 1.0;
        oldorient(NN, 8) = 1.0;

        int Nb = pots->nb;



        (*pots).nb = Nb + 1; //set the new pot

        obj->setdat(newdat);
        obj->setmom(oldmom);
        obj->setorientation(oldorient);
        obj->setangularmomenta(oldamom);

    }
    else{
        error("add 4 2 choice must be either 0 or 1");
    }

    }*/

    void NanotubeAssembly::setviscosity(double a)
    {
        double hdradius = 0.5;
        obj->setgamma(6. * pi * a * hdradius);
        obj->setgammar(8. * pi * a * hdradius * hdradius * hdradius);
    }

    void NanotubeAssembly::setkT(double a)
    {
        obj->setkT(a);
    }

    void NanotubeAssembly::setpots(ComboPatch &a)
    {
        //ComboPatch *q = a.clone();
        //delete pots;


        pots = (a).clone();

        // cout << (*(*pots).p)[0] << endl;
        // pausel();
        //pots =  a.clone();
        //(*pots).CreateFiles();
    }

    matrix<double> NanotubeAssembly::calculate_covariance(int Ns) {

        matrix<double> m(3,3);

        for(int i = 0 ; i < 3 ; i++) {
            for(int j = i ; j < 3 ; j++) {
                vector1<double> m1(Ns);
                vector1<double> m2(Ns);

                for(int k = 0 ; k < Ns ; k++) {
                    m1[k] = obj->getcoordinate(k,i);

                    m2[k] = obj->getcoordinate(k,j);
                }

                double mean1 = meanish(m1);
                double mean2 = meanish(m2);
                double tot = 0.0;
                for(int k = 0 ; k < Ns ; k++) {
                tot+=(m1[k] - mean1) * (m2[k] - mean2);
                }
                tot /= ((double)Ns-1);
            
                m(i,j) = tot;

                if(i != j) m(j,i) = tot; 
            }
        }
        return m;
    }

    void NanotubeAssembly::run(int runtime, int every, string strbase = "")
    {
        
        int ccc;

        int tf = ceil((double)runtime / (double)every);
        int number_of_digits = 0;
        do
        {
            ++number_of_digits;
            tf /= 10;
        } while (tf);

        matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

        matrix<int> *pairs = obj->calculatepairs(boxes, 3.5);

        WCAPotential wsa(3.0, 1.0, 0.0);

        int NN = obj->getN();

        cout << NN << endl;

        matrix<double> F(NN, 3);
        matrix<double> T(NN, 3);
        matrix<double> RT(NN, 6);
        matrix<double> zeromatrix(NN, 3);

        F = obj->calculateforces(*pairs, wsa);

        //cout << "ok to here" << endl;

        F += obj->calculateforces_external(conf);

        //cout << "trying to calculate this" << endl;
        obj->calculate_forces_and_torques3D(*pairs, *pots, F, T);

        //double coru = 1.0;
        // double nxtemp = 0.95;
        // double nytemp = 0.31225;
        // double nztemp = 0.0;
        // KernFrenkelOnePatch2 testpot(nxtemp, nytemp, -nztemp, nxtemp, nytemp, nztemp, 100., 2., pi / 3., 0.75);
        // obj->calculate_forces_and_torques3D(*pairs, testpot, F, T);

        //obj->create_random_forces(RT, RR);
        generate_uniform_random_matrix(RT);
        obj->create_forces_and_torques_sphere(F, T, RT);

        vector1<double> tottemp(6);

        for (int i = 0; i < runtime; i++)
        {
            //cout << i << endl;
            vector1<double> meas(6);
            // obj->measured_temperature(meas);
            // tottemp += meas;
            // cout << tottemp / (double)(i + 1) << endl;

            // cout << i << endl;
            if (i > 0 && i % 20 == 0)
            {
                // cout << "pairs recalculated" << endl;
                delete pairs;
                pairs = obj->calculatepairs(boxes, 3.5);
            }

            obj->advancemom_halfstep(F, T);

            obj->advance_pos();
            obj->rotate();

            F = obj->calculateforces(*pairs, wsa);

            F += obj->calculateforces_external(conf);
            //cout << obj->calculateforces_external(conf) << endl;
            //pausel();
            T.reset(0.0);

            obj->calculate_forces_and_torques3D(*pairs, *pots, F, T);

            // stringstream aa;
            // aa << setw(number_of_digits+1) << setfill('0') << (i / 1);
            // outfunc(T,"Tl_i="+aa.str());
            // outfunc(F, "Fl_i=" + aa.str());
            // obj->calculate_forces_and_torques3D(*pairs, *pots->potential_bundle[0], F, T);

            // obj->create_random_forces(RT, RR);
            generate_uniform_random_matrix(RT);
            obj->create_forces_and_torques_sphere(F, T, RT);

            // outfunc(T, "Tb_i=" + aa.str());
            // outfunc(F, "Fb_i=" + aa.str());

            obj->advancemom_halfstep(F, T);
            if (i % every == 0)
            {

                cout << i << endl;

                stringstream ss;

                ss << setw(number_of_digits) << setfill('0') << (i / every);

                matrix<double> orient = obj->getorientation();
                matrix<double> pos = obj->getdat();

                string poss = "pos";
                poss = poss + strbase;
                string oris = "orientation";
                oris = oris + strbase;

                vector<patchint> pairsbound = obj->calculate_bound_pairs(*pairs,*pots);

                poss += "_i=";
                oris += "_i=";

                string extension = ".csv";

                poss += ss.str();
                oris += ss.str();

                poss += extension;
                oris += extension;

                ofstream myfile;
                myfile.open(poss.c_str());

                ofstream myfile2;
                myfile2.open(oris.c_str());

                myfile <<= pos;
                //myfile2 <<= orient;
                for(int j = 0 ; j < pairsbound.size() ; j++) {
                myfile2 << pairsbound[j].particle_index1 << "," << pairsbound[j].particle_index2 << endl;
                }


                myfile.close();
                myfile2.close();

                //pausel();
            }
        }
    }

    void NanotubeAssembly::run_bending_modulus(int runtime, int every, double rate, string strbase = "")
    {

        int ccc;

        int tf = ceil((double)runtime / (double)every);
        int number_of_digits = 0;
        do
        {
            ++number_of_digits;
            tf /= 10;
        } while (tf);


        double ll = 50.;
        int total_beads = 2000;
        vector1<bool> pb(3, false);
        cube geo(ll+3., pb, 3);

        obj->setgeometry(geo);

        matrix<double> dat(total_beads, 3);

        double mp = 10.;

      

        for(int i  = 0 ; i < total_beads ; i++) {
            double rr1 = (double)rand() / (double)(RAND_MAX);
            double rr2 = (double)rand() / (double)(RAND_MAX);
            double rr3 = (double)rand() / (double)(RAND_MAX);
            dat(i,0) = 1+   ll*rr1;
            dat(i, 1) = 1 + ll * rr2;
            dat(i, 2) = 1 + ll * rr3;

            bool overlap=false;
            for(int j = 0 ; j < i ; j++ ) {
                
                overlap = geo.distance_less_than(dat, i, j, 1.2);
                
                if(overlap) {
                i--;
                break;
                }
            }
        }



        matrix<double> orientation(total_beads,9);
        for(int i = 0 ; i < total_beads ; i++) {
            orientation(i, 0)=1.;
            orientation(i, 4) = 1.;
            orientation(i, 8) = 1.;
        }
        matrix<double> mom(total_beads,3);
        matrix<double> angmom(total_beads,9);
       

        obj->setdat(dat);
        obj->setorientation(orientation);
        obj->setmom(mom);
        obj->setangularmomenta(angmom);

        vector1<double> lls(3);
        lls[0] = ll+3;
        lls[1] = ll+3;
        lls[2] = ll+3;

        planar_confinement conf2(lls,100.);

        int num2 = floor(ll / 4.);

        matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num2, ccc);

        matrix<int> *pairs = obj->calculatepairs(boxes, 3.5);

        WCAPotential wsa(3.0, 1.0, 0.0);

        int NN = obj->getN();

        cout << NN << endl;

  

        matrix<double> F(NN, 3);
        matrix<double> T(NN, 3);
        matrix<double> RT(NN, 6);
        matrix<double> zeromatrix(NN, 3);



        F = obj->calculateforces(*pairs, wsa);
   
        // double spring_constant = 100.;

        // F(total_beads - 1, 0) += -spring_constant * (obj->getcoordinate(total_beads - 1, 0) - mp);
        // F(total_beads - 1, 1) += - spring_constant*(obj->getcoordinate(total_beads-1, 1) - mp);
        // F(total_beads - 1, 2) += spring_constant * (ll + 1. - obj->getcoordinate(total_beads-1, 2));
        // F(0, 2) += -spring_constant * (obj->getcoordinate(0, 2));

        // F(0, 0) += -spring_constant * (obj->getcoordinate(0, 0) - mp);
        // F(0, 1) += -spring_constant * (obj->getcoordinate(0, 1) - mp);
        // cout << "ok to here" << endl;

        F += obj->calculateforces_external(conf2);


        // cout << "trying to calculate this" << endl;
        // cout << F << endl;
        
        obj->calculate_forces_and_torques3D(*pairs, *pots, F, T);
        // cout << F << endl;
        
        // double coru = 1.0;
        //  double nxtemp = 0.95;
        //  double nytemp = 0.31225;
        //  double nztemp = 0.0;
        //  KernFrenkelOnePatch2 testpot(nxtemp, nytemp, -nztemp, nxtemp, nytemp, nztemp, 100., 2., pi / 3., 0.75);
        //  obj->calculate_forces_and_torques3D(*pairs, testpot, F, T);

        // obj->create_random_forces(RT, RR);
        generate_uniform_random_matrix(RT);
        obj->create_forces_and_torques_sphere(F, T, RT);
        

        vector1<double> tottemp(6);

        for (int i = 0; i < runtime; i++)
        {
            ll -= rate;
            // cout << i << endl;
            vector1<double> meas(6);
            // obj->measured_temperature(meas);
            // tottemp += meas;
            // cout << tottemp / (double)(i + 1) << endl;

            // cout << i << endl;
            if (i > 0 && i % 20 == 0)
            {
                // cout << "pairs recalculated" << endl;
                delete pairs;
                pairs = obj->calculatepairs(boxes, 3.5);
            }
 
            obj->advancemom_halfstep(F, T);
            obj->advance_pos();
            obj->rotate();

            F = obj->calculateforces(*pairs, wsa);
            // F(total_beads - 1, 0) += -spring_constant * (obj->getcoordinate(total_beads - 1, 0) - mp);
            // F(total_beads - 1, 1) += -spring_constant * (obj->getcoordinate(total_beads - 1, 1) - mp);
            // F(total_beads - 1, 2) += spring_constant * (ll +1 - obj->getcoordinate(total_beads - 1, 2));
            // F(0, 2) += -spring_constant * (obj->getcoordinate(0, 2));
            // F(0, 0) += -spring_constant * (obj->getcoordinate(0, 0) - mp);
            // F(0, 1) += -spring_constant * (obj->getcoordinate(0, 1) - mp);

            lls[2] = ll + 1;
            conf2.setl(lls);
            F += obj->calculateforces_external(conf2);
            // cout << obj->calculateforces_external(conf) << endl;
            // pausel();
            T.reset(0.0);

            obj->calculate_forces_and_torques3D(*pairs, *pots, F, T);

            // stringstream aa;
            // aa << setw(number_of_digits+1) << setfill('0') << (i / 1);
            // outfunc(T,"Tl_i="+aa.str());
            // outfunc(F, "Fl_i=" + aa.str());
            // obj->calculate_forces_and_torques3D(*pairs, *pots->potential_bundle[0], F, T);

            // obj->create_random_forces(RT, RR);
            generate_uniform_random_matrix(RT);
            obj->create_forces_and_torques_sphere(F, T, RT);

            // outfunc(T, "Tb_i=" + aa.str());
            // outfunc(F, "Fb_i=" + aa.str());

            obj->advancemom_halfstep(F, T);
            if (i % every == 0)
            {

                cout << i << endl;

                stringstream ss;

                ss << setw(number_of_digits) << setfill('0') << (i / every);

                matrix<double> orient = obj->getorientation();
                matrix<double> pos = obj->getdat();

                string poss = "pos";
                poss = poss + strbase;
                string oris = "orientation";
                oris = oris + strbase;

                vector<patchint> pairsbound = obj->calculate_bound_pairs(*pairs, *pots);

                poss += "_i=";
                oris += "_i=";

                string extension = ".csv";

                poss += ss.str();
                oris += ss.str();

                poss += extension;
                oris += extension;

                ofstream myfile;
                myfile.open(poss.c_str());

                ofstream myfile2;
                myfile2.open(oris.c_str());

                myfile <<= pos;
                // myfile2 <<= orient;
                for (int j = 0; j < pairsbound.size(); j++)
                {
                    myfile2 << pairsbound[j].particle_index1 << "," << pairsbound[j].particle_index2 << endl;
                }

                myfile.close();
                myfile2.close();

                // pausel();
            }
        }
    }

    // void NanotubeAssembly::run_anneal(int runtime, int every, int cd, string strbase = "")
    // {
    //    // matrix<double> par = (*pots).params2;

    //     int nt2 = (*pots).nt;
    //     int nb2 = (*pots).nb;
    //     int nb_count = 0;
    //     // TetrahedralWithBivalent ax(par, nt2+nb2 - nb_count, nb_count);

    //     // this->setpots(ax);



    //     int ccc;

    //     int tf = ceil((double)runtime / (double)every);
    //     int number_of_digits = 0;
    //     do
    //     {
    //         ++number_of_digits;
    //         tf /= 10;
    //     } while (tf);

    //     matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    //     matrix<int> *pairs = obj->calculatepairs(boxes, 3.5);

    //     WCAPotential wsa(3.0, 1.0, 0.0);

    //     int NN = obj->getN();


    //     matrix<double> F(NN, 3);
    //     matrix<double> T(NN, 3);
    //     matrix<double> RT(NN, 6);
    //     matrix<double> zeromatrix(NN, 3);

    //     F = obj->calculateforces(*pairs, wsa);

    //     // cout << "ok to here" << endl;

    //     F += obj->calculateforces_external(conf);

    //     // cout << "trying to calculate this" << endl;
    //     obj->calculate_forces_and_torques3D(*pairs, *pots, F, T);

    //     // double coru = 1.0;
    //     //  double nxtemp = 0.95;
    //     //  double nytemp = 0.31225;
    //     //  double nztemp = 0.0;
    //     //  KernFrenkelOnePatch2 testpot(nxtemp, nytemp, -nztemp, nxtemp, nytemp, nztemp, 100., 2., pi / 3., 0.75);
    //     //  obj->calculate_forces_and_torques3D(*pairs, testpot, F, T);

    //     // obj->create_random_forces(RT, RR);
    //     generate_uniform_random_matrix(RT);
    //     obj->create_forces_and_torques_sphere(F, T, RT);

    //     vector1<double> tottemp(6);

    //     for (int i = 0; i < runtime; i++)
    //     {
    //         // cout << i << endl;
    //         vector1<double> meas(6);
    //         // obj->measured_temperature(meas);
    //         // tottemp += meas;
    //         // cout << tottemp / (double)(i + 1) << endl;

    //         // cout << i << endl;
    //         if (i > 0 && i % 20 == 0)
    //         {
    //             // cout << "pairs recalculated" << endl;
    //             delete pairs;
    //             pairs = obj->calculatepairs(boxes, 3.5);
    //         }

    //         obj->advancemom_halfstep(F, T);

    //         obj->advance_pos();
    //         obj->rotate();

    //         F = obj->calculateforces(*pairs, wsa);

    //         F += obj->calculateforces_external(conf);
    //         // cout << obj->calculateforces_external(conf) << endl;
    //         // pausel();
    //         T.reset(0.0);

    //         obj->calculate_forces_and_torques3D(*pairs, *pots, F, T);

    //         // stringstream aa;
    //         // aa << setw(number_of_digits+1) << setfill('0') << (i / 1);
    //         // outfunc(T,"Tl_i="+aa.str());
    //         // outfunc(F, "Fl_i=" + aa.str());
    //         // obj->calculate_forces_and_torques3D(*pairs, *pots->potential_bundle[0], F, T);

    //         // obj->create_random_forces(RT, RR);
    //         generate_uniform_random_matrix(RT);
    //         obj->create_forces_and_torques_sphere(F, T, RT);

    //         // outfunc(T, "Tb_i=" + aa.str());
    //         // outfunc(F, "Fb_i=" + aa.str());

    //         obj->advancemom_halfstep(F, T);
    //         if(i % cd == 0) {
    //             nb_count++;
    //             pots->change_nt(nt2 + nb2 - nb_count);

    //         }

    //         if (i % every == 0)
    //         {

    //             cout << i << endl;

    //             stringstream ss;

    //             ss << setw(number_of_digits) << setfill('0') << (i / every);

    //             matrix<double> orient = obj->getorientation();
    //             matrix<double> pos = obj->getdat();

    //             string poss = "pos";
    //             poss = poss + strbase;
    //             string oris = "orientation";
    //             oris = oris + strbase;

    //             poss += "_i=";
    //             oris += "_i=";

    //             string extension = ".csv";

    //             poss += ss.str();
    //             oris += ss.str();

    //             poss += extension;
    //             oris += extension;

    //             ofstream myfile;
    //             myfile.open(poss.c_str());

    //             ofstream myfile2;
    //             myfile2.open(oris.c_str());

    //             myfile <<= pos;
    //             myfile2 <<= orient;

    //             myfile.close();
    //             myfile2.close();

    //             // pausel();
    //         }
    //     }
    // }

    void ShellProperties::DoAnMC(double ll, bool FLAG = true)
    {

        int N = posi.getnrows();
        vector1<bool> vec(3);
        cube geo(ll, vec, 3);

        WCAPotential wsa(10.0, 1.0, 0.0);
        HarmonicPotential spr(k, rm);

        matrix<int> bindones(N, 7);
        vector1<int> iterators(N);

        for (int i = 0; i < par.getnrows(); i++)
        {
            int p1 = par(i, 0);
            int p2 = par(i, 1);

            int iter1 = iterators[p1];
            int iter2 = iterators[p2];

            bindones(p1, iter1) = p2;
            bindones(p2, iter2) = p1;

            iterators[p1]++;
            iterators[p2]++;
        }

        matrix<double> energy_per_bond(N, 7);
        for (int i = 0; i < N; i++)
        {
            int p1 = i;
            for (int j = 0; j < iterators[p1]; j++)
            {
                int p2 = bindones(i, j);

                double dis = geo.distance(posi, p1, p2);

                double en = wsa.energy(dis) + spr.energy(dis);

                energy_per_bond(i, j) = en;
            }
        }

        int No_Moves = 1000000;
        double step_size = 0.1;

        if (FLAG)
        {
            for (int i = 0; i < No_Moves; i++)
            {

                int choice = rand() % N; //choose a random particle

                double x = posi(choice, 0);
                double y = posi(choice, 1);
                double z = posi(choice, 2);

                double r1 = (double)rand() / (double)RAND_MAX;
                double r2 = (double)rand() / (double)RAND_MAX;
                double r3 = (double)rand() / (double)RAND_MAX;

                double newx = x + step_size * (2. * r1 - 1.);
                double newy = y + step_size * (2. * r2 - 1.);
                double newz = z + step_size * (2. * r3 - 1.);

                vector1<double> newpos(3);
                newpos[0] = newx;
                newpos[1] = newy;
                newpos[2] = newz;

                double befE = 0;
                for (int j = 0; j < iterators[choice]; j++)
                {
                    befE += energy_per_bond(choice, j);
                }

                double aftE = 0;
                vector1<double> es(iterators[choice]);
                for (int j = 0; j < iterators[choice]; j++)
                {
                    int p2 = bindones(choice, j);

                    double dis = geo.distance(newpos, posi.getrowvector(p2));

                    double en = wsa.energy(dis) + spr.energy(dis);

                    es[j] = en;
                    aftE += en;
                    //energy_per_bond(i, j) = en;
                }

                double rf = (double)(rand()) / (double)(RAND_MAX);
                double dE = aftE - befE;

                double kT = 1.;
                if (dE < 0)
                {
                    posi(choice, 0) = newx;
                    posi(choice, 1) = newy;
                    posi(choice, 2) = newz;
                    for (int j = 0; j < iterators[choice]; j++)
                    {
                        energy_per_bond(choice, j) = es[j];
                    }
                }
                if (exp(-dE / kT) > rf)
                {
                    posi(choice, 0) = newx;
                    posi(choice, 1) = newy;
                    posi(choice, 2) = newz;

                    for (int j = 0; j < iterators[choice]; j++)
                    {
                        energy_per_bond(choice, j) = es[j];
                    }
                }
                else
                {
                }

                // double gettheta = atan2(sqrt(x ^ 2 + y ^ 2), z);
                // double getphi = atan2(y,x);

                // double newgettheta = gettheta + r1;
                // double newgetphi = gettheta + r2;
            }
        }
        else{
        for (int i = 0; i < No_Moves; i++)
        {

            int choice = rand() % N; //choose a random particle
            double x = posi(choice, 0);
            double y = posi(choice, 1);
            double z = posi(choice, 2);

            double r1 = (double)rand() / (double)RAND_MAX;
            double r2 = (double)rand() / (double)RAND_MAX;
            // double r3 = (double)rand() / (double)RAND_MAX;
            double theta = atan2(sqrt(SQR(x)+SQR(y)),z);
            double phi = atan2(y,x);

            double r = sqrt(SQR(x)+SQR(y)+SQR(z));
            theta += step_size * (2. * r1 - 1.);
            phi += step_size * (2. * r2 - 1.);

            double newx = r * cos(phi) * sin(theta);
            double newy = r * sin(phi) * sin(theta);
            double newz = r * cos(theta);

            vector1<double> newpos(3);
            newpos[0] = newx;
            newpos[1] = newy;
            newpos[2] = newz;

            double befE = 0;
            for (int j = 0; j < iterators[choice]; j++)
            {
                befE += energy_per_bond(choice, j);
            }

            double aftE = 0;
            vector1<double> es(iterators[choice]);
            for (int j = 0; j < iterators[choice]; j++)
            {
                int p2 = bindones(choice, j);

                double dis = geo.distance(newpos, posi.getrowvector(p2));

                double en = wsa.energy(dis) + spr.energy(dis);

                es[j] = en;
                aftE += en;
                //energy_per_bond(i, j) = en;
            }

            double rf = (double)(rand()) / (double)(RAND_MAX);
            double dE = aftE - befE;

            double kT = 1.;
            if (dE < 0)
            {
                posi(choice, 0) = newx;
                posi(choice, 1) = newy;
                posi(choice, 2) = newz;
                for (int j = 0; j < iterators[choice]; j++)
                {
                    energy_per_bond(choice, j) = es[j];
                }
            }
            if (exp(-dE / kT) > rf)
            {
                posi(choice, 0) = newx;
                posi(choice, 1) = newy;
                posi(choice, 2) = newz;

                for (int j = 0; j < iterators[choice]; j++)
                {
                    energy_per_bond(choice, j) = es[j];
                }
            }
            else
            {
            }
            // double newx = x + step_size * (2. * r1 - 1.);
            // double newy = y + step_size * (2. * r2 - 1.);
            // double newz = z + step_size * (2. * r3 - 1.);
        

        }
    
    }
}
// struct ShellProperties {
//     matrix<int> par;
//     matrix<double> posi;
//     double k;
//     double rm;
// };

void NanotubeAssembly::run_with_real_surface(int runtime, int every, ShellProperties &myshell, matrix<double> &constantF, string strbase = "", bool cont = false)
{
    //constant F is a constant force applied on the system
    //WARNING, WE ARE NOT CHECKING WHETHER THE SHELL WE ARE ADDING IS NOT OVERLAPPING WITH THE ORIGINAL SYSTEM,
    //THE USER IS RESPONSIBLE FOR THIS
    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> bindingpairs =  myshell.par;

    matrix<int> quads =  myshell.quad;

    matrix<double> diss = myshell.bindingdis;

    int totnp = myshell.posi.getnrows();

 
    

    //Hard Sphere Forces
    WCAPotential wsa(3.0, 1.0, 0.0);
    HarmonicPotential spr(myshell.k,myshell.rm);

    // DoAnMC(ShellProperties);


    //Combine Our Data Into One
    int NN = obj->getN();
    if(cont) NN -= totnp;
    matrix<double> dat = obj->getdat();


    double ll = obj->getgeo().l;



    matrix<double> newdat(totnp+NN,3);
    vector1<int> p1(NN);
    if(cont) {
        newdat=dat;
        for (int i = totnp; i < totnp + NN; i++)
        {

                p1[i - totnp] = i;
            
        }
    }
    else{
        for(int i  = 0  ; i < totnp + NN ; i++) {
            if(i < totnp) {
                for(int j = 0  ; j < 3 ; j++)
                    newdat(i,j) = myshell.posi(i, j)+ll/2.;
            }
            else{
                for (int j = 0; j < 3; j++)
                    newdat(i, j) = dat(i-totnp, j);
                    p1[i-totnp] =  i;
            }
        }
    }




    
    obj->initialize(newdat);


    NN =  obj->getN(); //set new N
    

    matrix<int> *pairs = obj->calculatepairs_parallel(boxes, 3.5);
    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, p1, 3.5);


    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);



    F += obj->calculateforcesharmonic(bindingpairs, diss, myshell.k);



    F = obj->calculateforces(*pairs, wsa);

    F += constantF;
    // double val;
    // F.maxima(val);
    // cout << val << endl;
    // pausel();



    F += obj->calculateforcesdelauny(quads, myshell.kappa);

    // cout <<  obj->calculateforcesdelauny(quads, myshell.kappa) << endl;
    // cout << myshell.kappa << endl;
    // pausel();


    //cout << "ok to here" << endl;

    // F += obj->calculateforces_external(conf);

    //cout << "trying to calculate this" << endl;

    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T);


    //double coru = 1.0;
    // double nxtemp = 0.95;
    // double nytemp = 0.31225;
    // double nztemp = 0.0;
    // KernFrenkelOnePatch2 testpot(nxtemp, nytemp, -nztemp, nxtemp, nytemp, nztemp, 100., 2., pi / 3., 0.75);
    // obj->calculate_forces_and_torques3D(*pairs, testpot, F, T);

    //obj->create_random_forces(RT, RR);
    generate_uniform_random_matrix(RT);
    matrix<double> F2= F;
    obj->create_forces_and_torques_sphere(F, T, RT);

 
    vector1<double> tottemp(6);


    for (int i = 0; i < runtime; i++)
    {
        cout << i << endl;
        vector1<double> meas(6);
        // obj->measured_temperature(meas);
        // tottemp += meas;
        // cout << tottemp / (double)(i + 1) << endl;

        // cout << i << endl;
        if (i > 0 && i % 20 == 0)
        {
            // cout << "pairs recalculated" << endl;
            delete pairs;
            delete pairs_onlyb;
            // pairs = obj->calculatepairs(boxes, 3.5);
            pairs = obj->calculatepairs_parallel(boxes, 3.5);
            pairs_onlyb = obj->calculatepairs_parallel(boxes, p1, 3.5);

        }

/*         matrix<double> R(NN, 3);
        for (int i1 = 0; i1 < NN; i1++)
        {
            for (int j = 0; j < 3; j++)
            {
                R(i1, j) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
            }
        }
        F = obj->calculateforces(*pairs, wsa);

        (*obj).advance_mom(F, R);

        (*obj).advance_pos(); */



        obj->advancemom_halfstep(F, T);

        obj->advance_pos();


        obj->rotate();

        F = obj->calculateforces(*pairs, wsa);
        F += obj->calculateforcesharmonic(bindingpairs, diss, myshell.k);
        F += obj->calculateforcesdelauny(quads, myshell.kappa);
        F += constantF;

        // F += obj->calculateforces_external(conf);
        //cout << obj->calculateforces_external(conf) << endl;
        //pausel();
        T.reset(0.0);

        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T);

        // stringstream aa;
        // aa << setw(number_of_digits+1) << setfill('0') << (i / 1);
        // outfunc(T,"Tl_i="+aa.str());
        // outfunc(F, "Fl_i=" + aa.str());
        // obj->calculate_forces_and_torques3D(*pairs, *pots->potential_bundle[0], F, T);

        // obj->create_random_forces(RT, RR);
        generate_uniform_random_matrix(RT);
        matrix<double> F2 =  F;

        obj->create_forces_and_torques_sphere(F, T, RT);


        // outfunc(T, "Tb_i=" + aa.str());
        // outfunc(F, "Fb_i=" + aa.str());


        obj->advancemom_halfstep(F, T);

        if (i % every == 0)
        {

            cout << i << endl;

            stringstream ss;

            ss << setw(number_of_digits) << setfill('0') << (i / every);

            matrix<double> orient = obj->getorientation();
            matrix<double> pos = obj->getdat();

            string poss = "pos";
            poss = poss + strbase;
            string oris = "orientation";
            oris = oris + strbase;

            poss += "_i=";
            oris += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();

            poss += extension;
            oris += extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            myfile <<= pos;
            myfile2 <<= orient;

            myfile.close();
            myfile2.close();

            //pausel();
        }
    }
}


void NanotubeAssembly::run_with_real_surface_add_particles(int runtime, int every, ShellProperties &myshell, double prod, WeiM &c1 , string strbase = "")
{

    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> bindingpairs = myshell.par;

    matrix<int> quads = myshell.quad;

    matrix<double> diss = myshell.bindingdis;

    int totnp = myshell.posi.getnrows(); // total particles that are not patchy

    // Hard Sphere Forces
    HSPotential wsa(3.0, 1.0);
    HarmonicPotential spr(myshell.k, myshell.rm);

    int Nx = obj->getN(); // this defines the NN total, 
    //i.e. every particle that is going to be added to the simulation, 
    //but which might not be included yet
    int NN = Nx + totnp; //every particle
    
    matrix<double> dat = obj->getdat(); // get the data, remember that a lot of these values will be initialized
    //to zero to begin with


    double ll = obj->getgeo().l;

    matrix<double> newdat(NN, 3);
    for (int i = 0; i < NN; i++)
    {
        if (i < totnp)
        {
            for (int j = 0; j < 3; j++)
                newdat(i, j) = myshell.posi(i, j) + ll / 2.;
        }
        else
        {
            for (int j = 0; j < 3; j++)
                newdat(i, j) = dat(i - totnp, j);
        }
    }

    obj->initialize(newdat);


    //programmatic
    vector<int> indices_everything; //these are the indices of all the particles on the shell
    vector<int> indices_shell;
    indices_everything.reserve(NN);
    vector<int> indices_patchy; // the index of all the patchy particles
    for(int i = 0  ; i < totnp ; i++) {
        indices_everything.push_back(i);
        indices_shell.push_back(i);
    }
    indices_patchy.reserve(Nx);

    vector<int> indices_to_add;
    indices_to_add.reserve(Nx);
    vector<double> indices_weights;
    indices_weights.reserve(Nx);
    for (int i = totnp; i < NN; i++)
    {
        indices_to_add.push_back(i);
        double ww=1.;
        if(i<c1.M) {
            ww = c1.weight;
        }
        indices_weights.push_back(ww);
    }



    particle_adder vv;
    vv.set_indices(indices_to_add);
    vv.set_weights(indices_weights);
    sphere_vol vol;
    vol.r = (1./10.)*ll;
    vol.ll = ll;
    vol.c1 = ll/2.;
    vol.c2 = ll/2.;
    vol.c3 = ll/2.;
    vv.set_volume(vol);
    vv.set_rate(prod);

    
    matrix<int> *pairs = new matrix<int>;
    pairs = obj->calculatepairs_parallel(boxes,indices_everything, 3.5); //for the hard sphere repulsion

    matrix<int> *pairs_onlyb = new matrix<int>;

    pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5); // for the patchy particle binding

    //start with zero patchy particles in the simulation
    //indices_patchy.push_back(0);

    //define the full storage matrices of all the particles, despite the fact they won't all be utilized
    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    F = obj->calculateforces(*pairs, wsa); //calculate the forces due to hard sphere forces

    F += obj->calculateforcesharmonic(bindingpairs, diss, myshell.k);
    F += obj->calculateforcesdelauny(quads, myshell.kappa);

    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); //calculate the forces involved due to patchy

    generate_uniform_random_matrix(RT, indices_patchy); //only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT,indices_patchy,false); //only create torques and forces for patchy particles  

    for (int i = 0; i < runtime; i++)
    {
    cout << i << endl;
    if (i > 0 && i % 20 == 0)
    {
            // cout << "pairs recalculated" << endl;
            delete pairs;
            delete pairs_onlyb;
            // pairs = obj->calculatepairs(boxes, 3.5);
            
            //ADD THE PARTICLE ADDITION METHOD HERE
            bool dd = false;
            vector1<double> v1(3);
            int fi;


           vector1<double> me1 = meanish(obj->getdat(),indices_shell);

           sphere_vol vol2;
           vol2.r = (1. / 10.) * ll;
           vol2.ll = ll;
           vol2.c1 = me1[0];
           vol2.c2 = me1[1];
           vol2.c3 = me1[2];
           vv.set_volume(vol2);

           vv.add_p(*obj, indices_everything, dd, v1, fi);
            
           if (dd)
           {

                cout << "added: " << fi << endl;
                
                indices_everything.push_back(fi);
                
                indices_patchy.push_back(fi);
                obj->set_particle(v1, fi);
            }

            
            pairs = obj->calculatepairs_parallel(boxes,indices_everything, 3.5);

            pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5);

        }

    obj->advancemom_halfstep(F, T, indices_everything);
    obj->advance_pos(indices_everything);
    obj->rotate(indices_patchy); //only update the patchy particles with the rotate algorithm

    F = obj->calculateforces(*pairs, wsa);        // calculate the forces due to hard sphere forces

    F += obj->calculateforcesharmonic(bindingpairs, diss, myshell.k);

    F += obj->calculateforcesdelauny(quads, myshell.kappa);

    T.reset(0.0);


    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy


    generate_uniform_random_matrix(RT); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT); // only create torques and forces for patchy particles

    obj->advancemom_halfstep(F, T, indices_everything);

    if (i % every == 0 && i>0)
    {

            //cout << i << endl;

            stringstream ss;

            ss << setw(number_of_digits) << setfill('0') << (i / every);

            matrix<double> orient = obj->getorientation();
            matrix<double> pos = obj->getdat();

            string poss = "pos";
            poss = poss + strbase;
            string oris = "div";
            oris = oris + strbase;

            string elli = "or";
            elli = elli + strbase;

            poss += "_i=";
            oris += "_i=";
            //elli += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();
            elli += ss.str();

            poss += extension;
            oris += extension;
            elli += extension;

            string orie = "orient";
            orie+= extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            ofstream myfile3;
            myfile3.open(elli.c_str(),  std::ios_base::app);

            ofstream myfile4;
            myfile4.open(orie.c_str());

            myfile <<= pos;
            for(int ik  = 0 ; ik < indices_everything.size() ; ik++)
            myfile2 << indices_everything[ik] << endl;

            

            matrix<double> eig = calculate_covariance(totnp);

            myfile3 <<= obj->getorientation();

            myfile4 <<= obj->getorientation();



            myfile.close();
            myfile2.close();
            myfile3.close();
            myfile4.close();

            // pausel();
    }

    }
    delete pairs;
    delete pairs_onlyb;
}

template <typename T>
std::vector<T> complement(const std::vector<T> &v1, const std::vector<T> &v2)
{
    std::vector<T> result;
    std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(result));
    return result;
}

void NanotubeAssembly::run_with_real_surface_add_particles_continue(int runtime, int every, int starting_num, ShellProperties &myshell, double prod, WeiM &c1, matrix<double> &olddat, matrix<double> &oldori, vector1<int> &oldint, string strbase = "")
{

    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
    ++number_of_digits;
    tf /= 10;
    } while (tf);

    // int number_of_digits = 2; // as the code continues to iterate, we shall set it 

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> bindingpairs = myshell.par;

    matrix<int> quads = myshell.quad;

    matrix<double> diss = myshell.bindingdis;

    int totnp = myshell.posi.getnrows(); // total particles that are not patchy

    // Hard Sphere Forces
    HSPotential wsa(3.0, 1.0);
    HarmonicPotential spr(myshell.k, myshell.rm);

    // int Nx = obj->getN(); // this defines the NN total,
    // i.e. every particle that is going to be added to the simulation,
    // but which might not be included yet
    int NN = olddat.getnrows(); // every particle

    int Nx = NN - totnp;

    // cout << Nx << endl;
    // cout << NN << endl;
    // pausel();

    matrix<double> dat = olddat; // get the data, remember that a lot of these values will be initialized
    // to zero to begin with


    double ll = obj->getgeo().l;

    matrix<double> newdat(NN, 3); 
    for (int i = 0; i < NN; i++)
    {
        for (int j = 0; j < 3; j++)
        newdat(i, j) = dat(i, j);
    }


    obj->initialize(newdat);

    obj->setorientation(oldori);



    // programmatic
    vector<int> indices_everything; // these are the indices of all the particles on the shell
    vector<int> indices_shell;
    indices_everything.reserve(NN);
    vector<int> indices_patchy; // the index of all the patchy particles

    indices_patchy.reserve(Nx);
    for (int i = 0; i < totnp; i++)
    {
    indices_everything.push_back(i);
    indices_shell.push_back(i);
    }
    vector<int> already_added; // the particles that have already been added
    for(int i = totnp ; i < oldint.getsize() ; i++) {
        already_added.push_back(oldint[i]);
        indices_everything.push_back(oldint[i]);
        indices_patchy.push_back(oldint[i]);
    }
    

    vector<int> indices_to_add_temp;
    indices_to_add_temp.reserve(Nx);
    vector<double> indices_weights;
    indices_weights.reserve(Nx);
    for (int i = totnp; i < NN; i++)
    {
    indices_to_add_temp.push_back(i);
    }

    // for(int i  = 0 ; i < already_added.size() ; i++){
    // int index_which = already_added[i];
    // remove_at(indices_weights, index_which);
    // remove_at(indices_to_add, index_which);
    // }


    std::sort(already_added.begin(),already_added.end());

    vector<int> indices_to_add  = complement(indices_to_add_temp,already_added);



   
   
    // set_difference(indices_to_add_temp.begin(),indices_to_add_temp.end(),already_added.begin(),already_added.end(),std::back_inserter(indices_to_add));

    for (int i = 0; i < indices_to_add.size(); i++)
    {
    //indices_to_add_temp.push_back(i);
    int index1 = indices_to_add[i];
    double ww = 1.;
    if (index1 < c1.M)
    {
        ww = c1.weight;
    }
    indices_weights.push_back(ww);
    }

    // for (size_t i = 0; i < already_added.size(); ++i)
    // {
    // std::cout << already_added[i] << ",";
    // }
    // cout << endl;
    // cout << endl;

    // for (size_t i = 0; i < indices_to_add.size(); ++i)
    // {
    // std::cout << indices_to_add[i] << ",";
    // }
    // cout << endl;

    // cout << already_added.size() << endl;
    // cout << indices_weights.size() << endl;
    // cout << indices_to_add.size() << endl;
    // pausel();

    particle_adder vv;
    vv.set_indices(indices_to_add);
    vv.set_weights(indices_weights);
    sphere_vol vol;
    vol.r = (1. / 10.) * ll;
    vol.ll = ll;
    vol.c1 = ll / 2.;
    vol.c2 = ll / 2.;
    vol.c3 = ll / 2.;
    vv.set_volume(vol);
    vv.set_rate(prod);

    matrix<int> *pairs = obj->calculatepairs_parallel(boxes, indices_everything, 3.5); // for the hard sphere repulsion

    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5); // for the patchy particle binding

    vector<patchint> pat = obj->calculate_patch_list(*pairs_onlyb, *pots);



    // start with zero patchy particles in the simulation
    // indices_patchy.push_back(0);

    // define the full storage matrices of all the particles, despite the fact they won't all be utilized
    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    F = obj->calculateforces(*pairs, wsa); // calculate the forces due to hard sphere forces

    F += obj->calculateforcesharmonic(bindingpairs, diss, myshell.k);
    F += obj->calculateforcesdelauny(quads, myshell.kappa);

    obj->calculate_forces_and_torques3D(pat, *pots, F, T); // calculate the forces involved due to patchy

    // obj->calculate_forces_and_torques3D(pat, *pots, F2, T2); // calculate the forces involved due to patchy



    generate_uniform_random_matrix(RT, indices_patchy); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT, indices_patchy, false); // only create torques and forces for patchy particles


    for (int i = 0; i < runtime; i++)
    {
    cout << i << endl;
    if (i > 0 && i % 20 == 0)
    {
            // cout << "pairs recalculated" << endl;
            delete pairs;
            delete pairs_onlyb;
            // pairs = obj->calculatepairs(boxes, 3.5);

            // ADD THE PARTICLE ADDITION METHOD HERE
            bool dd = false;
            vector1<double> v1(3);
            int fi;

            vector1<double> me1 = meanish(obj->getdat(), indices_shell);

            sphere_vol vol2;
            vol2.r = (1. / 10.) * ll;
            vol2.ll = ll;
            vol2.c1 = me1[0];
            vol2.c2 = me1[1];
            vol2.c3 = me1[2];
            vv.set_volume(vol2);

            vv.add_p(*obj, indices_everything, dd, v1, fi);

            if (dd)
            {
            cout << "added: " << fi << endl;

            indices_everything.push_back(fi);

            indices_patchy.push_back(fi);

            obj->set_particle(v1, fi);

            }

            pairs = obj->calculatepairs_parallel(boxes, indices_everything, 3.5);

            pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5);

            pat = obj->calculate_patch_list(*pairs_onlyb, *pots);
    }

    obj->advancemom_halfstep(F, T, indices_everything);
    obj->advance_pos(indices_everything);
    obj->rotate(indices_patchy); // only update the patchy particles with the rotate algorithm

    F = obj->calculateforces(*pairs, wsa); // calculate the forces due to hard sphere forces

    F += obj->calculateforcesharmonic(bindingpairs, diss, myshell.k);
    F += obj->calculateforcesdelauny(quads, myshell.kappa);

    T.reset(0.0);

    obj->calculate_forces_and_torques3D(pat, *pots, F, T); // calculate the forces involved due to patchy

    generate_uniform_random_matrix(RT); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT); // only create torques and forces for patchy particles

    obj->advancemom_halfstep(F, T, indices_everything);

    if (i>0 && i % every == 0)
    {

            //cout << i << endl;

            stringstream ss;

            ss << setw(number_of_digits) << setfill('0') << starting_num  +(i / every);
            //cout << "done 1" << endl;
            matrix<double> orient = obj->getorientation();
            matrix<double> pos = obj->getdat();

            string poss = "pos";
            poss = poss + strbase;
            string oris = "div";
            oris = oris + strbase;

            string elli = "eigs";
            elli = elli + strbase;

            poss += "_i=";
            oris += "_i=";
            // elli += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();
            // elli += ss.str();

            poss += extension;
            oris += extension;
            elli += extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            ofstream myfile3;
            myfile3.open(elli.c_str(), std::ios_base::app);

            myfile <<= pos;
            for (int ik = 0; ik < indices_everything.size(); ik++)
            myfile2 << indices_everything[ik] << endl;

            matrix<double> eig = calculate_covariance(totnp);

            //cout << "done 2" << endl;
            myfile3 <<= eig;

            myfile.close();
            myfile2.close();
            myfile3.close();

            string orie = "orient"; //we want to overwrite it, as it is a lot of data without much utility
            orie += extension;
            ofstream myfile4;
            myfile4.open(orie.c_str());
            myfile4 <<= obj->getorientation();
            //cout << "done 3" << endl;
            // pausel();
    }
    }

    delete pairs;
    delete pairs_onlyb;
}

    void NanotubeAssembly::run_with_real_surface_remove_add_particles_continue(int runtime, int every, int starting_num, ShellProperties &myshell, double prod, double k1, double k2, WeiM &c1, matrix<double> &olddat, matrix<double> &oldori, vector1<int> &oldint, string strbase = "")
    {

        int ccc;

        int tf = ceil((double)runtime / (double)every);
        int number_of_digits = 0;
        do
        {
            ++number_of_digits;
            tf /= 10;
        } while (tf);

        // int number_of_digits = 2; // as the code continues to iterate, we shall set it

        matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

        matrix<int> bindingpairs = myshell.par; // we can get all the binding pairs of the elastic shell

        int totnp = myshell.posi.getnrows(); // total particles that are not patchy

        // Hard Sphere Forces
        HSPotential wsa(3.0, 1.0);
        HarmonicPotential spr(myshell.k, myshell.rm);

        int NN = olddat.getnrows(); // every particle

        int Nx = NN - totnp;

        // cout << Nx << endl;
        // cout << NN << endl;
        // pausel();

        matrix<double> dat = olddat; // get the data, remember that a lot of these values will be initialized
        // to zero to begin with

        double ll = obj->getgeo().l;

        matrix<double> newdat(NN, 3);
        for (int i = 0; i < NN; i++)
        {
            for (int j = 0; j < 3; j++)
                newdat(i, j) = dat(i, j);
        }

        obj->initialize(newdat);

        obj->setorientation(oldori);

        cout << NN << endl;
        cout << totnp << endl;
        cout << Nx << endl;
        // programmatic
        vector<int> indices_everything; // these are the indices of all the particles on the shell
        vector<int> indices_shell;
        indices_everything.reserve(NN);
        vector<int> indices_patchy; // the index of all the patchy particles

        indices_patchy.reserve(Nx);
        for (int i = 0; i < totnp; i++)
        {
            indices_everything.push_back(i);
            indices_shell.push_back(i);
        }
        vector<int> already_added; // the particles that have already been added
        for (int i = totnp; i < oldint.getsize(); i++)
        {
            already_added.push_back(oldint[i]);
            indices_everything.push_back(oldint[i]);
            indices_patchy.push_back(oldint[i]);
        }

        //we want a scheme that removes things from indices_patchy and adds them back at some rate too
        vector<int> indices_removed;

        matrix<int> *pairs = obj->calculatepairs_parallel(boxes, indices_everything, 3.5); // for the hard sphere repulsion

        matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5); // for the patchy particle binding

        matrix<double> F(NN, 3);
        matrix<double> Fs(NN, 3);
        matrix<double> T(NN, 3);
        matrix<double> RT(NN, 6);
        matrix<double> zeromatrix(NN, 3);

        F = obj->calculateforces(*pairs, wsa); // calculate the forces due to hard sphere forces

        F += obj->calculateforces(bindingpairs, spr); // calculate the forces involved due to elastic shell

        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

        generate_uniform_random_matrix(RT, indices_patchy); // only generate random torques for the patchy particles

        obj->create_forces_and_torques_sphere(F, T, RT, indices_patchy, false); // only create torques and forces for patchy particles

        for (int i = 0; i < runtime; i++)
        {
            cout << i << endl;
            if (i > 0 && i % 20 == 0)
            {
                // cout << "pairs recalculated" << endl;
                delete pairs;
                delete pairs_onlyb;
                // pairs = obj->calculatepairs(boxes, 3.5);

                double total_reaction_rate=k1*indices_patchy.size()+k2*indices_removed.size();
                double time_between = 1./total_reaction_rate;
                double r1 = ((double)rand() / (double)(RAND_MAX));
                // ADD THE PARTICLE ADDITION METHOD HERE
                //either nothing, inactivate an active particle, or reactivate an inactivated particle
                cout << total_reaction_rate << endl;



                if(r1<total_reaction_rate) {
                    double forward=k1*indices_patchy.size()/total_reaction_rate;
                    double backward = 1-forward;

                    double r2 = ((double)rand() / (double)(RAND_MAX));
                    if(r2 < forward) {
                         //choose a random particle, remove it from list indices_patchy and add it to vector indicies_removed
                    
                        int randelement = rand() % indices_patchy.size();

                        auto it = indices_patchy.begin();
                        std::advance(it,randelement);
                        indices_removed.push_back(*it);
                        indices_patchy.erase(it);

                        // cout << *it << endl;
                        }

                    else{
                        int randelement = rand() % indices_removed.size();

                        auto it = indices_removed.begin();
                        std::advance(it, randelement);
                        indices_patchy.push_back(*it);
                        indices_removed.erase(it);

                    }
                }

                

                pairs = obj->calculatepairs_parallel(boxes, indices_everything, 3.5);

                pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5);
            }

            obj->advancemom_halfstep(F, T, indices_everything);
            obj->advance_pos(indices_everything);
            obj->rotate(indices_patchy); // only update the patchy particles with the rotate algorithm

            F = obj->calculateforces(*pairs, wsa); // calculate the forces due to hard sphere forces

            F += obj->calculateforces(bindingpairs, spr); // calculate the forces involved due to elastic shell

            T.reset(0.0);

            obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

            generate_uniform_random_matrix(RT); // only generate random torques for the patchy particles

            obj->create_forces_and_torques_sphere(F, T, RT); // only create torques and forces for patchy particles

            obj->advancemom_halfstep(F, T, indices_everything);

            if (i > 0 && i % every == 0)
            {

                cout << i << endl;

                stringstream ss;

                ss << setw(number_of_digits) << setfill('0') << starting_num + (i / every);
                cout << "done 1" << endl;
                matrix<double> orient = obj->getorientation();
                matrix<double> pos = obj->getdat();

                string poss = "pos";
                poss = poss + strbase;
                string oris = "div";
                oris = oris + strbase;

                string elli = "eigs";
                elli = elli + strbase;

                poss += "_i=";
                oris += "_i=";
                // elli += "_i=";

                string extension = ".csv";

                poss += ss.str();
                oris += ss.str();
                // elli += ss.str();

                poss += extension;
                oris += extension;
                elli += extension;

                ofstream myfile;
                myfile.open(poss.c_str());

                ofstream myfile2;
                myfile2.open(oris.c_str());

                ofstream myfile3;
                myfile3.open(elli.c_str(), std::ios_base::app);

                myfile <<= pos;
                for (int ik = 0; ik < indices_everything.size(); ik++)
                    myfile2 << indices_everything[ik] << endl;

                // matrix<double> eig = calculate_covariance(totnp);

                for (int ik = 0; ik < indices_patchy.size()-1; ik++)
                    myfile3 << indices_patchy[ik] << ",";

                myfile3 << indices_patchy[indices_patchy.size() - 1] << endl;
                if(indices_removed.size() >0) {
                for (int ik = 0; ik < indices_removed.size()-1; ik++)
                    myfile3 << indices_removed[ik] << ",";

                myfile3 << indices_removed[indices_removed.size() - 1] << endl;
                }
                myfile.close();
                myfile2.close();
                myfile3.close();

                string orie = "orient"; // we want to overwrite it, as it is a lot of data without much utility
                orie += extension;
                ofstream myfile4;
                myfile4.open(orie.c_str());
                myfile4 <<= obj->getorientation();
                cout << "done 3" << endl;
                // pausel();
            }
        }

        delete pairs;
        delete pairs_onlyb;
    }
 
 
    void NanotubeAssembly::run_add_particles(int runtime, int every, double prod, string strbase = "")
    {

        int ccc;

        int tf = ceil((double)runtime / (double)every);
        int number_of_digits = 0;
        do
        {
            ++number_of_digits;
            tf /= 10;
        } while (tf);

        double ll2 = obj->getgeo().l;
        vector1<bool> pb(3, false);
        cube geo(ll2, pb, 3);
        int num2 = floor(ll2 / 4.);
        obj->setgeometry(geo);

        matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num2, ccc);

        // matrix<int> bindingpairs = myshell.par; // we can get all the binding pairs of the elastic shell

        // int totnp = myshell.posi.getnrows(); // total particles that are not patchy

        // Hard Sphere Forces
        HSPotential wsa(3.0, 1.0);
        // HarmonicPotential spr(myshell.k, myshell.rm);

        int Nx = obj->getN(); // this defines the NN total,
        // i.e. every particle that is going to be added to the simulation,
        // but which might not be included yet
        int totnp = Nx;
        int NN = Nx; // + totnp; // every particle

        // programmatic
        vector<int> indices_everything; // these are the indices of all the particles on the shell

        indices_everything.reserve(NN);
        vector<int> indices_patchy; // the index of all the patchy particles
        for (int i = 0; i < totnp; i++)
        {
            indices_everything.push_back(i);
        }
        indices_patchy.reserve(Nx);

        vector<int> indices_to_add;
        indices_to_add.reserve(Nx);
        vector<double> indices_weights;
        indices_weights.reserve(Nx);
        for (int i = 0; i < NN; i++)
        {
            indices_to_add.push_back(i);
            double ww = 1.;

            indices_weights.push_back(ww);
        }

        particle_adder vv;
        vv.set_indices(indices_to_add);
        vv.set_weights(indices_weights);
        sphere_vol vol;

        vol.r = (1. / 10.) * ll;
        vol.ll = ll;
        vol.c1 = ll / 2.;
        vol.c2 = ll / 2.;
        vol.c3 = ll / 2.;
        vv.set_volume(vol);
        vv.set_rate(prod);

        // start with zero patchy particles in the simulation
        // indices_patchy.push_back(0);

        // define the full storage matrices of all the particles, despite the fact they won't all be utilized
        matrix<double> F(NN, 3);
        matrix<double> T(NN, 3);
        matrix<double> RT(NN, 6);
        matrix<double> zeromatrix(NN, 3);

        matrix<int> *pairs = new matrix<int>;
        pairs = obj->calculatepairs_parallel(boxes, indices_everything, 3.5); // for the hard sphere repulsion

        matrix<int> *pairs_onlyb = new matrix<int>;

        pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5); // for the patchy particle binding

        F = obj->calculateforces(*pairs, wsa); // calculate the forces due to hard sphere forces

        cout << F << endl;
        cout << NN << endl;
        pausel();
        // F += obj->calculateforces(bindingpairs, spr); // calculate the forces involved due to elastic shell
        F += obj->calculateforces_external(conf);

        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

        generate_uniform_random_matrix(RT, indices_patchy); // only generate random torques for the patchy particles

        obj->create_forces_and_torques_sphere(F, T, RT, indices_patchy, false); // only create torques and forces for patchy particles

        for (int i = 0; i < runtime; i++)
        {
            cout << i << endl;
            if (i > 0 && i % 20 == 0)
            {
                // cout << "pairs recalculated" << endl;
                delete pairs;
                delete pairs_onlyb;
                // pairs = obj->calculatepairs(boxes, 3.5);

                // ADD THE PARTICLE ADDITION METHOD HERE
                bool dd = false;
                vector1<double> v1(3);
                int fi;

                sphere_vol vol2;
                vol2.r = (1. / 10.) * ll;
                vol2.ll = ll;
                vol2.c1 = ll / 2;
                vol2.c2 = ll / 2;
                vol2.c3 = ll / 2;
                vv.set_volume(vol2);

                vv.add_p(dd, fi);

                if (dd)
                {
                    cout << "added: " << fi << endl;
                    // indices_everything.push_back(fi);
                    indices_patchy.push_back(fi);
                    // obj->set_particle(v1, fi);
                }

                pairs = obj->calculatepairs_parallel(boxes, indices_everything, 3.5);

                pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5);
            }

            obj->advancemom_halfstep(F, T, indices_everything);
            obj->advance_pos(indices_everything);
            obj->rotate(indices_patchy); // only update the patchy particles with the rotate algorithm

            F = obj->calculateforces(*pairs, wsa); // calculate the forces due to hard sphere forces

            // F += obj->calculateforces(bindingpairs, spr); // calculate the forces involved due to elastic shell
            F += obj->calculateforces_external(conf);
            T.reset(0.0);

            obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

            // if(pairs_onlyb->getnrows() > 0 ) {
            // cout << *pairs_onlyb << endl;
            // vector1<double> v1 = F.getrowvector(pairs_onlyb->gpcons(0, 0));
            // vector1<double> v2 = F.getrowvector(pairs_onlyb->gpcons(0, 1));

            // cout << v1 << endl;
            // cout << v2 << endl;
            // if(abs(v1[0])>1E-10) {
            //     pausel();
            // }

            // }

            generate_uniform_random_matrix(RT, indices_patchy); // only generate random torques for the patchy particles

            obj->create_forces_and_torques_sphere(F, T, RT, indices_patchy, false); // only create torques and forces for patchy particles

            obj->advancemom_halfstep(F, T, indices_everything);

            if (i % every == 0)
            {

                cout << i << endl;

                stringstream ss;

                ss << setw(number_of_digits) << setfill('0') << (i / every);

                matrix<double> orient = obj->getorientation();
                matrix<double> pos = obj->getdat();

                string poss = "pos";
                poss = poss + strbase;
                string oris = "div";
                oris = oris + strbase;

                poss += "_i=";
                oris += "_i=";

                string extension = ".csv";

                poss += ss.str();
                oris += ss.str();

                poss += extension;
                oris += extension;

                ofstream myfile;
                myfile.open(poss.c_str());

                ofstream myfile2;
                myfile2.open(oris.c_str());

                myfile <<= pos;
                for (int ik = 0; ik < indices_patchy.size(); ik++)
                    myfile2 << indices_patchy[ik] << endl;

                myfile.close();
                myfile2.close();

                // pausel();
            }
        }
        delete pairs;
        delete pairs_onlyb;
    }

#endif /* NANOTUBE_CPP */
