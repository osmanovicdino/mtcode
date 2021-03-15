#ifndef NANOTUBE_CPP
#define NANOTUBE_CPP

NanotubeAssembly::NanotubeAssembly(double rmax, int N) {


    double ll =  rmax;

    num = floor(ll/ 4.);

    obj = new LangevinNVTR;

    BivalentPatch temp_patch(100., 1.3, pi / 4.); //initialize abstract base class to be a one patch system

    //pots = new ComboPatch;

    this->setpots(temp_patch);



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

    double viscosity = 1.0;
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

            //cout << i << endl;

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
        }
    }
}

#endif /* NANOTUBE_CPP */
