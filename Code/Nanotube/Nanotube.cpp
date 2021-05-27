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


    double ll =  rmax;

    num = floor(ll/ 4.);

    obj = new LangevinNVTR;

    BivalentPatch temp_patch(100., 1.3, pi / 4.); //initialize abstract base class to be a one patch system

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
            bool cond3 = i >=12 && j >= 12;



            if( cond1  || cond2 ) {
                params(iter, 0) = 10.0;
                params(iter, 1) = 1.4;
                params(iter, 2) = 0.4;
                iter++;
            }
            else if(cond3) {
                params(iter, 0) = 30.0;
                params(iter, 1) = 1.2;
                params(iter, 2) = 0.4;
                iter++;
            }
            else{
                params(iter, 0) = 0.0;
                params(iter, 1) = 1.4;
                params(iter, 2) = 0.3;
                iter++;
            }
        }
    }

    outfunc(params,"lig");

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

    double dphi = 0.2;
    double dtheta =0.1;
    vector1<double> p1 = B(1, dtheta + pid / 2., -dphi);
    vector1<double> p2 = B(1, pid / 2., -dphi);
    vector1<double> p3 = B(1, -dtheta + pid / 2., -dphi);

    vector1<double> p4 = B(1, dtheta + pid / 2.,(2.*pid/3.)  + dphi);
    vector1<double> p5 = B(1, pid / 2., (2. * pid / 3.) + dphi );
    vector1<double> p6 = B(1, -dtheta + pid / 2., (2. * pid / 3.) + dphi );

    vector1<double> p7 = B(1, dtheta + pid / 2., dphi);
    vector1<double> p8 = B(1, pid / 2., dphi);
    vector1<double> p9 = B(1, -dtheta + pid / 2., dphi);

    vector1<double> p10 = B(1, dtheta + pid / 2., (2. * pid / 3.) - dphi);
    vector1<double> p11 = B(1, pid / 2., (2. * pid / 3.) - dphi);
    vector1<double> p12 = B(1, -dtheta + pid / 2., (2. * pid / 3.) - dphi);

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

    orient(0,0) = p1[0];
    orient(0,1) = p1[1];
    orient(0,2) = p1[2];

    orient(1,0) = p2[0];
    orient(1,1) = p2[1];
    orient(1,2) = p2[2];

    orient(2,0) = p3[0];
    orient(2,1) = p3[1];
    orient(2,2) = p3[2];

    orient(3,0) = p4[0];
    orient(3,1) = p4[1];
    orient(3,2) = p4[2];

    orient(4,0) = p5[0];
    orient(4,1) = p5[1];
    orient(4,2) = p5[2];

    orient(5,0) = p6[0];
    orient(5,1) = p6[1];
    orient(5,2) = p6[2];

    orient(6,0) = p7[0];
    orient(6,1) = p7[1];
    orient(6,2) = p7[2];

    orient(7,0) = p8[0];
    orient(7,1) = p8[1];
    orient(7,2) = p8[2];

    orient(8,0) = p9[0];
    orient(8,1) = p9[1];
    orient(8,2) = p9[2];

    orient(9,0) = p10[0];
    orient(9,1) = p10[1];
    orient(9,2) = p10[2];

    orient(10,0) = p11[0];
    orient(10,1) = p11[1];
    orient(10,2) = p11[2];

    orient(11,0) = p12[0];
    orient(11,1) = p12[1];
    orient(11,2) = p12[2];

    // orient(12,0) = 0.;
    // orient(12,1) = 0.;
    // orient(12,2) = 1.0;

    // orient(13, 0) = 0.;
    // orient(13, 1) = 0.;
    // orient(13, 2) = -1.0;

    GeneralPatch c(vec1, numb, params, orient, true);


    BivalentPatch c2(10.0, 1.4,0.5);

    this->setpots(c2);



    spherical_confinement_3D conf2(rmax/2.,1.0,ll/2.);
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

    for (int i = 0; i < pp; i++)
    {
        for (int j = 0; j < pp; j++)
        {
            for (int k = 0; k < pp; k++)
            {
                double x = 0.5 + i;
                double y = 0.5 + j;
                double z = 0.5 + k;
                if (SQR(x - ll / 2.) + SQR(y - ll / 2.) + SQR(z - ll / 2.) < SQR(0.9 * rmax/2.) )
                {
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

    LangevinNVTR b(geo);

    b.initialize(dat);

    double kT = 1.0;
    double dt = 0.005;
    b.setdt(dt);

    double viscosity = 1.;
    double hdradius = 0.5;

    b.setgamma(6. * pi * viscosity * hdradius);
    b.setgammar(8. * pi * viscosity * hdradius * hdradius * hdradius);

    b.setkT(kT);

    *obj = b;


}

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
    pots = (a).clone();


    // cout << (*(*pots).p)[0] << endl;
    // pausel();
    //pots =  a.clone();
    (*pots).CreateFiles();
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

#endif /* NANOTUBE_CPP */
