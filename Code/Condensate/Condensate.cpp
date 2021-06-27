

Condensate::Condensate(double ll, int N)  {
    num = floor(ll/4.);
    ls = ll;

    obj = new LangevinNVTR;

    SingPatch temp_patch(100.,2.,pi/3.); //initialize abstract base class to be a one patch system

    this->setpots(temp_patch);

    BindingModelSingle temp_bind(1.0,0.0);

    this->setBindingModel(temp_bind);

    vector1<bool> pb(3,true);
    cube geo(ll,pb,3);

    matrix<double> dat(N,3);


    int pp =  ceil(cbrt((double)N)) +1;
    double binsize = (ll)/(double)pp;




    vector<double> possible_pos_x;
    vector<double> possible_pos_y;
    vector<double> possible_pos_z;

    for(int i = 0 ; i < pp ; i++) {
        for(int j = 0  ; j < pp ; j++) {
            for(int k = 0 ; k < pp ; k++) {
                double x = 0.5*binsize + i*binsize;
                double y = 0.5*binsize + j*binsize;
                double z = 0.5*binsize + k*binsize;
                possible_pos_x.push_back(x);
                possible_pos_y.push_back(y);
                possible_pos_z.push_back(z);
            }
        }
    }



   

    for(int i = 0 ; i < N ; i++) {
        int randint  =  rand() % (possible_pos_x.size());
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
    

    b.setgamma(6.*pi*viscosity*hdradius);
    b.setgammar(8.*pi*viscosity*hdradius*hdradius*hdradius);

    b.setkT(kT);

    *obj = b;

}

void Condensate::setviscosity(double a)
{
    double hdradius = 0.5;
    obj->setgamma(6. * pi * a * hdradius);
    obj->setgammar(8. * pi * a * hdradius * hdradius * hdradius);
}

void Condensate::setkT(double a) {
    obj->setkT(a);
}

void Condensate::setpots(ComboPatch &a) {
    ComboPatch *q = a.clone();
    pots = q;
    (*pots).CreateFiles();
}

void Condensate::setBindingModel(AbstractBindingModel &a) { 
    AbstractBindingModel *q = a.clone();
    bm = q;
}

void Condensate::run(int runtime, int every, string strbase = "")
{
    int ccc;





    int tf = ceil( (double)runtime / (double)every);
    int number_of_digits = 0;
    do {++number_of_digits;
    tf /= 10;
    } while(tf);



    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *pairs = obj->calculatepairs(boxes, 3.5);

    WCAPotential wsa(3.0,1.0,0.0);

    int NN = obj->getN();

    matrix<double> F(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN,3);


    F = obj->calculateforces(*pairs, wsa);
    
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
        obj->measured_temperature(meas);
        tottemp += meas;
        cout << tottemp/(double)(i+1) << endl;


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

void Condensate::run_singlebond(int runtime, int every, string strbase = "")
{

    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);



    int NN = obj->getN();


    BinaryBindStore bbs;

    
    int nh  = (*pots).get_total_patches(NN);

 

    vector1<bool>isbound(nh);

    vector1<int> boundto(nh);

    bbs.isbound = isbound;
    bbs.boundto =  boundto;

    num = floor(ls / (2. * size_mol));
    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);



    matrix<int> *pairs = obj->calculatepairs(boxes, 3.5*size_mol);


    WCAPotential wsa(1.0, size_mol, 0.0);



    matrix<double> F(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);



    obj->calculateforces(*pairs, wsa);

    vector<patchint> opairs;
    opairs.reserve(NN * 10 * 10);
    vector<int> runs_diff;
    runs_diff.reserve(NN * 20);
    opairs = obj->calculate_patch_list(*pairs, *pots);
    runs_diff = adjacency(opairs);

    //cout << "fi" << endl;
    //cout << "fi" << endl;
    obj->calculate_forces_and_torques3D_onlyone_nonlets(opairs, runs_diff, *pots, bbs , *bm, F, T);
    //cout << "fi2" << endl;


    generate_uniform_random_matrix(RT);
    obj->create_forces_and_torques_sphere(F, T, RT);
    //vector1<double> tottemp(6);

    for (int i = 0; i < runtime; i++)
    {

        cout << i << endl;
        // vector1<double> meas(6);
        // obj->measured_temperature(meas);
        // tottemp += meas;
        // cout << tottemp / (double)(i + 1) << endl;
        if (i > 0 && i % 20 == 0)
        {
           // cout << "pairs recalculated" << endl;
            delete pairs;
            pairs = obj->calculatepairs(boxes, 3.5);
        }
        if(i > 0 && i %10 ==0 ){
            opairs = obj->calculate_patch_list(*pairs, *pots);
            runs_diff = adjacency(opairs);
        }

       
        //obj->measured_temperature();

      
        obj->advancemom_halfstep(F, T);

        obj->advance_pos();

        obj->rotate();

        F = obj->calculateforces(*pairs, wsa);

        T.reset(0.0);

        obj->calculate_forces_and_torques3D_onlyone_nonlets(opairs, runs_diff, *pots, bbs, *bm, F, T);

        generate_uniform_random_matrix(RT);
        obj->create_forces_and_torques_sphere(F, T, RT);

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
            string bins =  "bindings";
            bins =  bins + strbase;

            poss += "_i=";
            oris += "_i=";
            bins += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();
            bins += ss.str();

            poss += extension;
            oris += extension;
            bins += extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            ofstream myfile3;
            myfile3.open(bins.c_str());

            myfile <<= pos;
            myfile2 << setprecision(10) <<= orient;
            myfile3 <<= bbs.isbound;
            myfile3 << "\n";
            myfile3 <<= bbs.boundto;


            myfile.close();
            myfile2.close();
            myfile3.close();
        }
    }
}

void Condensate::run_singlebond_different_sizes(int runtime, int every, int div, double size1, double size2, string strbase = "")
{

    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    int NN = obj->getN();

    vector1<int> p1(div);
    vector1<int> p2(NN-div);
    for(int i = 0  ; i < div ; i++) {
        p1[i] = i;
    }
    for(int i = div ; i < NN ; i++) {
        p2[i-div] = i;
    }

    BinaryBindStore bbs;

    int nh = (*pots).get_total_patches(NN);

    vector1<bool> isbound(nh);

    vector1<int> boundto(nh);

    bbs.isbound = isbound;
    bbs.boundto = boundto;

    
    num = floor(ls / (3.*size1) );

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *pairsp1 = obj->calculatepairs(boxes, p1, 3.5*size1);

    matrix<int> *pairsp2 = obj->calculatepairs(boxes, p2, 3.5*size2);

    matrix<int> *pairsp1p2 = obj->calculatepairs(boxes, p1, p2, 3.5*(size1+size2)/2.);

    WCAPotential wsa1(3.0, size1, 0.0);
    WCAPotential wsa2(3.0, (size1+size2)/2., 0.0);
    WCAPotential wsa3(3.0, size2, 0.0);

    matrix<double> F(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    F = obj->calculateforces(*pairsp1, wsa1);
    F += obj->calculateforces(*pairsp2, wsa3);
    F += obj->calculateforces(*pairsp1p2, wsa2);

    /*make combined pairs*/
    int npp = (pairsp1->getNsafe())+(pairsp2->getNsafe())+(pairsp1p2->getNsafe());
    matrix<int> pairs;
    matrix<int> pairstemp(npp,2);
    int iter = 0;
    for(int j = 0 ; j < pairsp1->getNsafe() ; j++) {
        pairstemp(iter, 0) = pairsp1->operator()(j, 0);
        pairstemp(iter, 1) = pairsp1->operator()(j, 1);
        iter++;
    }
    for (int j = 0; j < pairsp2->getNsafe(); j++)
    {
        pairstemp(iter, 0) = pairsp2->operator()(j, 0);
        pairstemp(iter, 1) = pairsp2->operator()(j, 1);
        iter++;
    }
    for (int j = 0; j < pairsp1p2->getNsafe(); j++)
    {
        pairstemp(iter, 0) = pairsp1p2->operator()(j, 0);
        pairstemp(iter, 1) = pairsp1p2->operator()(j, 1);
        iter++;
    }

    pairs = pairstemp;
    vector<patchint> opairs;
    opairs.reserve(NN * 10 * 10);
    vector<int> runs_diff;
    runs_diff.reserve(NN * 20);
    opairs = obj->calculate_patch_list(pairs, *pots);
    runs_diff = adjacency(opairs);

    obj->calculate_forces_and_torques3D_onlyone_nonlets(opairs, runs_diff, *pots, bbs, *bm, F, T);
    generate_uniform_random_matrix(RT);
    obj->create_forces_and_torques_sphere(F, T, RT);
    //vector1<double> tottemp(6);

    for (int i = 0; i < runtime; i++)
    {

        cout << i << endl;
        //obj->measured_temperature();

        // vector1<double> meas(6);
        // obj->measured_temperature(meas);
        // tottemp += meas;
        // cout << tottemp / (double)(i + 1) << endl;
        if (i > 0 && i % 20 == 0)
        {
            // cout << "pairs recalculated" << endl;
            delete pairsp1;
            delete pairsp2;
            delete pairsp1p2;

            pairsp1 = obj->calculatepairs(boxes, p1, 3.5 * size1);

            pairsp2 = obj->calculatepairs(boxes, p2, 3.5 * size2);

            pairsp1p2 = obj->calculatepairs(boxes, p1, p2, 3.5 * (size1 + size2) / 2.);

            npp = (pairsp1->getNsafe()) + (pairsp2->getNsafe()) + (pairsp1p2->getNsafe());
            
            matrix<int> pairstemp2(npp, 2);
            int iter = 0;
            for (int j = 0; j < pairsp1->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp1->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp1->operator()(j, 1);
                iter++;
            }
            for (int j = 0; j < pairsp2->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp2->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp2->operator()(j, 1);
                iter++;
            }
            for (int j = 0; j < pairsp1p2->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp1p2->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp1p2->operator()(j, 1);
                iter++;
            }

            pairs = pairstemp2;
        }

        if(i>0 && i %10 ==0 ) {
            opairs = obj->calculate_patch_list(pairs, *pots);
            runs_diff = adjacency(opairs);
        }

        //obj->measured_temperature();

        obj->advancemom_halfstep(F, T);


        obj->advance_pos();



        obj->rotate();

        F = obj->calculateforces(*pairsp1, wsa1);
        F += obj->calculateforces(*pairsp2, wsa3);
        F += obj->calculateforces(*pairsp1p2, wsa2);

        T.reset(0.0);

        obj->calculate_forces_and_torques3D_onlyone_nonlets(opairs, runs_diff,  *pots, bbs, *bm, F, T);
        generate_uniform_random_matrix(RT);
        obj->create_forces_and_torques_sphere(F, T, RT);

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
            string bins = "bindings";
            bins = bins + strbase;

            poss += "_i=";
            oris += "_i=";
            bins += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();
            bins += ss.str();

            poss += extension;
            oris += extension;
            bins += extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            ofstream myfile3;
            myfile3.open(bins.c_str());

            myfile <<= pos;
            myfile2 << setprecision(10) <<= orient;
            myfile3 <<= bbs.isbound;
            myfile3 << "\n";
            myfile3 <<= bbs.boundto;

            myfile.close();
            myfile2.close();
            myfile3.close();
        }
    }
}

void Condensate::run_singlebond_different_sizes_continue(int runtime, int every, int div, double size1, double size2, int starval, BinaryBindStore &bbs2, string strbase = "")
{

    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    int NN = obj->getN();

   
    vector1<int> p1(div);
    vector1<int> p2(NN - div);
    for (int i = 0; i < div; i++)
    {
        p1[i] = i;
    }
    for (int i = div; i < NN; i++)
    {
        p2[i - div] = i;
    }

    BinaryBindStore bbs = bbs2;

    int nh = (*pots).get_total_patches(NN);

    if(nh != bbs2.boundto.getsize() ) error("passing bindings incorrectly in run_singlebond_different_sizes_continue");
    // vector1<bool> isbound(nh);

    // vector1<int> boundto(nh);

    // bbs.isbound = isbound;
    // bbs.boundto = boundto;

    num = floor(ls / (3. * size1));

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *pairsp1 = obj->calculatepairs(boxes, p1, 3.5 * size1);

    matrix<int> *pairsp2 = obj->calculatepairs(boxes, p2, 3.5 * size2);

    matrix<int> *pairsp1p2 = obj->calculatepairs(boxes, p1, p2, 3.5 * (size1 + size2) / 2.);

    WCAPotential wsa1(3.0, size1, 0.0);
    WCAPotential wsa2(3.0, (size1+size2)/2., 0.0);
    WCAPotential wsa3(3.0, size2, 0.0);

    matrix<double> F(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    F = obj->calculateforces(*pairsp1, wsa1);

    F += obj->calculateforces(*pairsp2, wsa3);

    F += obj->calculateforces(*pairsp1p2, wsa2);



    /*make combined pairs*/
    int npp = (pairsp1->getNsafe()) + (pairsp2->getNsafe()) + (pairsp1p2->getNsafe());
    matrix<int> pairs;
    matrix<int> pairstemp(npp, 2);
    int iter = 0;
    for (int j = 0; j < pairsp1->getNsafe(); j++)
    {
        pairstemp(iter, 0) = pairsp1->operator()(j, 0);
        pairstemp(iter, 1) = pairsp1->operator()(j, 1);
        iter++;
    }
    for (int j = 0; j < pairsp2->getNsafe(); j++)
    {
        pairstemp(iter, 0) = pairsp2->operator()(j, 0);
        pairstemp(iter, 1) = pairsp2->operator()(j, 1);
        iter++;
    }
    for (int j = 0; j < pairsp1p2->getNsafe(); j++)
    {
        pairstemp(iter, 0) = pairsp1p2->operator()(j, 0);
        pairstemp(iter, 1) = pairsp1p2->operator()(j, 1);
        iter++;
    }

    pairs = pairstemp;
    vector<patchint> opairs;
    opairs.reserve(NN * 10 * 10);
    vector<int> runs_diff;
    runs_diff.reserve(NN * 20);
    opairs = obj->calculate_patch_list(pairs, *pots);
    runs_diff = adjacency(opairs);

    obj->calculate_forces_and_torques3D_onlyone_nonlets(opairs, runs_diff ,*pots, bbs, *bm, F, T);
    generate_uniform_random_matrix(RT);
    obj->create_forces_and_torques_sphere(F, T, RT);
    //vector1<double> tottemp(6);


    for (int i = 0; i < runtime; i++)
    {

        cout << i << endl;

        // vector1<double> meas(6);
        // obj->measured_temperature(meas);
        // tottemp += meas;
        // cout << tottemp / (double)(i + 1) << endl;
        if (i > 0 && i % 10 == 0)
        {
            // cout << "pairs recalculated" << endl;


            delete pairsp1;
            delete pairsp2;
            delete pairsp1p2;

            pairsp1 = obj->calculatepairs(boxes, p1, 3.5 * size1);

            pairsp2 = obj->calculatepairs(boxes, p2, 3.5 * size2);

            pairsp1p2 = obj->calculatepairs(boxes, p1, p2, 3.5 * (size1 + size2) / 2.);

            npp = (pairsp1->getNsafe()) + (pairsp2->getNsafe()) + (pairsp1p2->getNsafe());

            matrix<int> pairstemp2(npp, 2);
            int iter = 0;
            for (int j = 0; j < pairsp1->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp1->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp1->operator()(j, 1);
                iter++;
            }
            for (int j = 0; j < pairsp2->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp2->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp2->operator()(j, 1);
                iter++;
            }
            for (int j = 0; j < pairsp1p2->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp1p2->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp1p2->operator()(j, 1);
                iter++;
            }

            pairs = pairstemp2;
        }
        if(i>0 && i%10==0) {
            opairs = obj->calculate_patch_list(pairs, *pots);
            runs_diff = adjacency(opairs);
        }

        //obj->measured_temperature();

        obj->advancemom_halfstep(F, T);

        obj->advance_pos();

        //matrix<double> tempor = obj->getorientation();

        

        obj->rotate();


        // tempor =  absmatrix(tempor - obj->getorientation());
        // double valb;
        // tempor.maxima(valb);
        // cout << valb << endl;
        // pausel();

        F = obj->calculateforces(*pairsp1, wsa1);
        F += obj->calculateforces(*pairsp2, wsa3);
        F += obj->calculateforces(*pairsp1p2, wsa2);

        T.reset(0.0);

        obj->calculate_forces_and_torques3D_onlyone_nonlets(opairs, runs_diff, *pots, bbs, *bm, F, T);
        generate_uniform_random_matrix(RT);
        obj->create_forces_and_torques_sphere(F, T, RT);

        obj->advancemom_halfstep(F, T);

        if (i % every == 0)
        {

            //cout << i << endl;

            stringstream ss;

            ss << setw(number_of_digits) << setfill('0') << starval + (i / every);

            cout << "outputting file" << endl;
            cout << ss.str() << endl;

           // pausel();


            matrix<double> orient = obj->getorientation();
            matrix<double> pos = obj->getdat();

            string poss = "pos";
            poss = poss + strbase;
            string oris = "orientation";
            oris = oris + strbase;
            string bins = "bindings";
            bins = bins + strbase;

            poss += "_i=";
            oris += "_i=";
            bins += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();
            bins += ss.str();

            poss += extension;
            oris += extension;
            bins += extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            ofstream myfile3;
            myfile3.open(bins.c_str());

            myfile <<= pos;
            myfile2 << setprecision(10) <<= orient;
            myfile3 <<= bbs.isbound;
            myfile3 << "\n";
            myfile3 <<= bbs.boundto;

            myfile.close();
            myfile2.close();
            myfile3.close();
        }
    }
}

void Condensate::run_singlebond_different_sizes_continue_thetalist(int runtime, int every, int div, double size1, double size2, int starval, BinaryBindStore &bbs2, string strbase = "")
{

    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    int NN = obj->getN();

    vector1<int> p1(div);
    vector1<int> p2(NN - div);
    for (int i = 0; i < div; i++)
    {
        p1[i] = i;
    }
    for (int i = div; i < NN; i++)
    {
        p2[i - div] = i;
    }

    BinaryBindStore bbs = bbs2;

    int nh = (*pots).get_total_patches(NN);

    if (nh != bbs2.boundto.getsize())
        error("passing bindings incorrectly in run_singlebond_different_sizes_continue");
    // vector1<bool> isbound(nh);

    // vector1<int> boundto(nh);

    // bbs.isbound = isbound;
    // bbs.boundto = boundto;

    num = floor(ls / (3. * size1));

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *pairsp1 = obj->calculatepairs(boxes, p1, 3.5 * size1);

    matrix<int> *pairsp2 = obj->calculatepairs(boxes, p2, 3.5 * size2);

    matrix<int> *pairsp1p2 = obj->calculatepairs(boxes, p1, p2, 3.5 * (size1 + size2) / 2.);

    WCAPotential wsa1(3.0, size1, 0.0);
    WCAPotential wsa2(3.0, (size1 + size2) / 2., 0.0);
    WCAPotential wsa3(3.0, size2, 0.0);

    matrix<double> F(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    F = obj->calculateforces(*pairsp1, wsa1);

    F += obj->calculateforces(*pairsp2, wsa3);

    F += obj->calculateforces(*pairsp1p2, wsa2);

    /*make combined pairs*/
    int npp = (pairsp1->getNsafe()) + (pairsp2->getNsafe()) + (pairsp1p2->getNsafe());
    matrix<int> pairs;
    matrix<int> pairstemp(npp, 2);
    int iter = 0;
    for (int j = 0; j < pairsp1->getNsafe(); j++)
    {
        pairstemp(iter, 0) = pairsp1->operator()(j, 0);
        pairstemp(iter, 1) = pairsp1->operator()(j, 1);
        iter++;
    }
    for (int j = 0; j < pairsp2->getNsafe(); j++)
    {
        pairstemp(iter, 0) = pairsp2->operator()(j, 0);
        pairstemp(iter, 1) = pairsp2->operator()(j, 1);
        iter++;
    }
    for (int j = 0; j < pairsp1p2->getNsafe(); j++)
    {
        pairstemp(iter, 0) = pairsp1p2->operator()(j, 0);
        pairstemp(iter, 1) = pairsp1p2->operator()(j, 1);
        iter++;
    }

    pairs = pairstemp;



    obj->calculate_forces_and_torques3D_onlyone_nonlets(pairs, *pots, bbs, *bm, F, T);

    generate_uniform_random_matrix(RT);
    obj->create_forces_and_torques_sphere(F, T, RT);
    //vector1<double> tottemp(6);

    for (int i = 0; i < runtime; i++)
    {

        cout << i << endl;

        // vector1<double> meas(6);
        // obj->measured_temperature(meas);
        // tottemp += meas;
        // cout << tottemp / (double)(i + 1) << endl;
        if (i > 0 && i % 20 == 0)
        {
            // cout << "pairs recalculated" << endl;

            delete pairsp1;
            delete pairsp2;
            delete pairsp1p2;

            pairsp1 = obj->calculatepairs(boxes, p1, 3.5 * size1);

            pairsp2 = obj->calculatepairs(boxes, p2, 3.5 * size2);

            pairsp1p2 = obj->calculatepairs(boxes, p1, p2, 3.5 * (size1 + size2) / 2.);

            npp = (pairsp1->getNsafe()) + (pairsp2->getNsafe()) + (pairsp1p2->getNsafe());

            matrix<int> pairstemp2(npp, 2);
            int iter = 0;
            for (int j = 0; j < pairsp1->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp1->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp1->operator()(j, 1);
                iter++;
            }
            for (int j = 0; j < pairsp2->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp2->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp2->operator()(j, 1);
                iter++;
            }
            for (int j = 0; j < pairsp1p2->getNsafe(); j++)
            {
                pairstemp2(iter, 0) = pairsp1p2->operator()(j, 0);
                pairstemp2(iter, 1) = pairsp1p2->operator()(j, 1);
                iter++;
            }

            pairs = pairstemp2;
        }


        //obj->measured_temperature();

        obj->advancemom_halfstep(F, T);

        obj->advance_pos();

        //matrix<double> tempor = obj->getorientation();

        obj->rotate();

        // tempor =  absmatrix(tempor - obj->getorientation());
        // double valb;
        // tempor.maxima(valb);
        // cout << valb << endl;
        // pausel();

        F = obj->calculateforces(*pairsp1, wsa1);
        F += obj->calculateforces(*pairsp2, wsa3);
        F += obj->calculateforces(*pairsp1p2, wsa2);

        T.reset(0.0);
        
        obj->calculate_forces_and_torques3D_onlyone_nonlets(pairs, *pots, bbs, *bm, F, T);
        //obj->calculate_forces_and_torques3D_onlyone_nonlets(opairs,runs_diff, *pots, bbs3, *bm, F, T);
       
     
        generate_uniform_random_matrix(RT);
        obj->create_forces_and_torques_sphere(F, T, RT);

        obj->advancemom_halfstep(F, T);

        if (i % every == 0)
        {

            //cout << i << endl;

            stringstream ss;

            ss << setw(number_of_digits) << setfill('0') << starval + (i / every);

            cout << "outputting file" << endl;
            cout << ss.str() << endl;

            // pausel();

            matrix<double> orient = obj->getorientation();
            matrix<double> pos = obj->getdat();

            string poss = "pos";
            poss = poss + strbase;
            string oris = "orientation";
            oris = oris + strbase;
            string bins = "bindings";
            bins = bins + strbase;

            poss += "_i=";
            oris += "_i=";
            bins += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();
            bins += ss.str();

            poss += extension;
            oris += extension;
            bins += extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            ofstream myfile3;
            myfile3.open(bins.c_str());

            myfile <<= pos;
            myfile2 << setprecision(10) <<= orient;
            myfile3 <<= bbs.isbound;
            myfile3 << "\n";
            myfile3 <<= bbs.boundto;

            myfile.close();
            myfile2.close();
            myfile3.close();
        }
    }
}