{

    double r1 = ((double)rand() / (double)(RAND_MAX));
    double r2 = ((double)rand() / (double)(RAND_MAX));
    double r3 = ((double)rand() / (double)(RAND_MAX));
    double x1 = ll / 2. + (0.9 * myrmax / 2.) * (2. * r1 - 1.);
    double y1 = ll / 2. + (0.9 * myrmax / 2.) * (2. * r2 - 1.);
    double z1 = ll / 2. + (0.9 * myrmax / 2.) * (2. * r3 - 1.); 

    //choose a random position

    if (which == 0)
    {
        // add a 4

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
            returned double y2 = olddat(i, 1);
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

        // swap the rows

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

        (*pots).nt = Nt + 1; // set the new pot

        // for (int i = 2048; i < NN + 1; i++)
        //     cout << newdat.getrowvector(i) << endl;
        // cout << endl;

        // pausel();

        obj->setdat(newdat);
        obj->setmom(oldmom);
        obj->setorientation(oldorient);
        obj->setangularmomenta(oldamom);
    }
    else if (which == 1)
    {
        // add a 2

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

        (*pots).nb = Nb + 1; // set the new pot

        obj->setdat(newdat);
        obj->setmom(oldmom);
        obj->setorientation(oldorient);
        obj->setangularmomenta(oldamom);
    }
    else
    {
        error("add 4 2 choice must be either 0 or 1");
    }
}

void NanotubeAssembly::run_with_real_surface_add_particles(int runtime, int every, ShellProperties &myshell, double prod, string strbase = "")
{

    // WARNING, WE ARE NOT CHECKING WHETHER THE SHELL WE ARE ADDING IS NOT OVERLAPPING WITH THE ORIGINAL SYSTEM,
    // THE USER IS RESPONSIBLE FOR THIS
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

    int totnp = myshell.posi.getnrows(); // total particles that are not patchy

    // Hard Sphere Forces
    HSPotential wsa(3.0, 1.0);
    HarmonicPotential spr(myshell.k, myshell.rm);

    // DoAnMC(ShellProperties);

    // Combine Our Data Into One
    int NN = obj->getN();
    matrix<double> dat = obj->getdat();

    double ll = obj->getgeo().l;
    matrix<double> newdat(totnp + NN, 3);
    for (int i = 0; i < totnp + NN; i++)
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

    NN = obj->getN();            // set new N
    vector1<int> p1(NN - totnp); // binders only

    for (int ik = totnp; ik < NN; ik++)
    {
        p1[ik - totnp] = ik;
    }
    matrix<int> *pairs = obj->calculatepairs_parallel(boxes, 3.5);
    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, p1, 3.5);

    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    F = obj->calculateforces(*pairs, wsa);

    // double val;
    // F.maxima(val);
    // cout << val << endl;
    // pausel();

    F += obj->calculateforces(bindingpairs, spr);

    // cout << "ok to here" << endl;

    // F += obj->calculateforces_external(conf);

    // cout << "trying to calculate this" << endl;

    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T);

    // double coru = 1.0;
    //  double nxtemp = 0.95;
    //  double nytemp = 0.31225;
    //  double nztemp = 0.0;
    //  KernFrenkelOnePatch2 testpot(nxtemp, nytemp, -nztemp, nxtemp, nytemp, nztemp, 100., 2., pi / 3., 0.75);
    //  obj->calculate_forces_and_torques3D(*pairs, testpot, F, T);

    // obj->create_random_forces(RT, RR);
    generate_uniform_random_matrix(RT);
    matrix<double> F2 = F;

    obj->create_forces_and_torques_sphere(F, T, RT, vector1<int>(), true);

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
            double r1 = (double)rand() / (double)RAND_MAX;

            if (r1 < prod)
            {
                int ad = rand() % 10;
                int ad2;
                if (ad < 9)
                    ad2 = 1;
                else
                    ad2 = 0;
                this->add_particle42(ad2);
                int Ng = obj->getN();
                F.resize(Ng, 3);
                Fs.resize(Ng, 3);
                T.resize(Ng, 3);
                RT.resize(Ng, 6);
                zeromatrix.resize(Ng, 3);
            }

            NN = obj->getN(); // set new N
            vector1<int> p2(NN - totnp);
            for (int ik = totnp; ik < NN; ik++)
            {
                p2[ik - totnp] = ik;
            }

            pairs = obj->calculatepairs_parallel(boxes, 3.5);
            pairs_onlyb = obj->calculatepairs_parallel(boxes, p2, 3.5);
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

        F += obj->calculateforces(bindingpairs, spr);
        // F += obj->calculateforces_external(conf);
        // cout << obj->calculateforces_external(conf) << endl;
        // pausel();
        T.reset(0.0);

        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T);

        // stringstream aa;
        // aa << setw(number_of_digits+1) << setfill('0') << (i / 1);
        // outfunc(T,"Tl_i="+aa.str());
        // outfunc(F, "Fl_i=" + aa.str());
        // obj->calculate_forces_and_torques3D(*pairs, *pots->potential_bundle[0], F, T);

        // obj->create_random_forces(RT, RR);
        generate_uniform_random_matrix(RT);

        matrix<double> F2 = F;

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
            myfile2 << pots->nt;

            myfile.close();
            myfile2.close();

            // pausel();
        }
    }
}