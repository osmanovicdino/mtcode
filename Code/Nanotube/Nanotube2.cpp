#ifndef NANOTUBE2_CPP
#define NANOTUBE2_CPP


bool test_overlap(matrix<double> &a, int v, cube &geo, double overlap_distance) {

    // this tests the overlaps of polymers which are of size 3;
    bool overlap = false;
    if(v == 0 ) return overlap;

    // v is the index of the current polymer
    for(int k = 0 ; k < v ; k++) {
        for(int j = 0 ; j < 3 ; j++) {
            for(int i = 0  ; i < 3 ; i++) {
                double dis = geo.distance(a.getrowvector(k*3+j),a.getrowvector(v*3+i));
                if(dis < overlap_distance) {
                    overlap = true;
                    goto overfound; //skip the rest of the check
                }
            }
        }


    }
    overfound:
    return overlap;
}

void NanotubeAssembly::run_no_patchy(int runtime, int every, ShellProperties &myshell, double prod, double ks, double ks2, string strbase = "")
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

    //placeholder
    HarmonicPotential pol(ks,0.0);

    BendingPotential bp(ks2,0.);



    int Nx = obj->getN(); // this defines the NN total,
    // i.e. every particle that is going to be added to the simulation,
    // but which might not be included yet
    int NN = Nx + totnp; // every particle

    matrix<double> dat = obj->getdat(); // get the data, remember that a lot of these values will be initialized
    // to zero to begin with

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


    

    //create the polymers
    int NpolB=40; //number of polymers
   
    matrix<double> origpospolymers(NpolB*3,3);

    for(int k = 0 ; k < NpolB ; k++) {
        

        double r1 = (double)rand()/(double)RAND_MAX;
        double r2 = (double)rand() / (double)RAND_MAX;
        double r3 = (double)rand() / (double)RAND_MAX;

        double theta = acos(2*r1-1);
        double phi =  2*pi*r2;
        double psi = 2*pi*r3;


        matrix<double> rotmatrix = rotationmatrix(theta,phi,psi);

        vector1<double> p1(3);
        vector1<double> p2(3);
        vector1<double> p3(3);
        p2[0]=1.1;
        p3[0]=2.2;

        vector1<double> v1(3,ll/2.);

        double rmax =12.; //max r to displace;


        // the polymer is defined to be a line
        double r11 = (double)rand() / (double)RAND_MAX;
        double r12 = (double)rand() / (double)RAND_MAX;
        double r13 = (double)rand() / (double)RAND_MAX;

        double rm2 = rmax*r13;
        double theta2 = acos(2 * r11 - 1);
        double phi2 = 2 * pi * r12;
        // apply a translation and a random rotation to the polymer, as long as it is inside the shell
        //ROTATE
        p1 = rotmatrix * (p1);
        p2 = rotmatrix * (p2);
        p3 = rotmatrix * (p3);

        v1[0] += rm2 * cos(phi2) * sin(theta2); // rmax * (-1.+2. * (double)rand() / (double)RAND_MAX);
        v1[1] += rm2 * sin(phi2) * sin(theta2); // rmax * (-1.+2. * (double)rand() / (double)RAND_MAX);
        v1[2] += rm2 * cos(theta2);             // rmax * (-1.+2. * (double)rand() / (double)RAND_MAX);

        // TRANSLATE
        p1 += v1;
        p2 += v1;
        p3 += v1;



        for(int l = 0 ; l < 3 ; l++) {
        origpospolymers(k * 3 + 0, l)=p1[l];
        origpospolymers(k * 3 + 1, l) = p2[l];
        origpospolymers(k * 3 + 2, l) = p3[l];
        }
        cube testgeo = obj->getgeo();
         bool overlaptest = test_overlap(origpospolymers, k, testgeo, 1.1);



        if(overlaptest) k--; //if an overlap is found, try again

       


    }


    
    //Add these polymers to the system

    
    for (int i = totnp; i < totnp+origpospolymers.getnrows(); i++)
    {
        for (int j = 0; j < 3; j++)
            newdat(i, j) = origpospolymers(i - totnp, j);
        
    }

    obj->initialize(newdat);


    // programmatic
    vector<int> indices_everything; // these are the indices of all the particles on the shell
    vector<int> indices_shell;
    indices_everything.reserve(NN);
    vector<int> indices_patchy; // the index of all the patchy particles
    for (int i = 0; i < totnp; i++)
    {
        indices_everything.push_back(i);
        indices_shell.push_back(i);
    }
    for(int i = totnp ; i < totnp +NpolB*3 ; i++) {
        indices_everything.push_back(i);
    }
    indices_patchy.reserve(Nx);


    Nx = Nx - NpolB*3;
    vector<int> indices_to_add;
    indices_to_add.reserve(Nx);
    vector<double> indices_weights;
    indices_weights.reserve(Nx);
    for (int i = totnp+NpolB*3; i < NN; i++)
    {
        indices_to_add.push_back(i);

    }

    // particle_adder vv;
    // vv.set_indices(indices_to_add);
    // vv.set_weights(indices_weights);
    // sphere_vol vol;
    // vol.r = (1. / 10.) * ll;
    // vol.ll = ll;
    // vol.c1 = ll / 2.;
    // vol.c2 = ll / 2.;
    // vol.c3 = ll / 2.;
    // vv.set_volume(vol);
    // vv.set_rate(prod);

    matrix<int> *pairs = new matrix<int>;
    pairs = obj->calculatepairs_parallel(boxes, indices_everything, 3.5); // for the hard sphere repulsion

    vector<mdtriplet> triplets;
    vector<mdpair> polpairs;

    vector<vector<int> > allpolymers(NpolB);

    for(int i =  0 ; i < NpolB ; i++) {


        for(int j = 0  ; j < 3 ; j++)
        allpolymers[i].push_back(totnp+ i*3+j);

        mdtriplet mdt(totnp + i * 3, totnp + i * 3 + 1, totnp + i * 3 + 2);
        mdpair mdp1(totnp + i * 3, totnp + i * 3 + 1);
        mdpair mdp2(totnp + i * 3 + 1, totnp + i * 3 + 2);
        triplets.push_back(mdt);
        polpairs.push_back(mdp1);
        polpairs.push_back(mdp2);
    }


    // matrix<int> *pairs_onlyb = new matrix<int>;

    // pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5); // for the patchy particle binding

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
    F += obj->calculateforces(polpairs,pol);
    F += obj->calculateforces_threebody(triplets, bp);

    // obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

    generate_uniform_random_matrix(RT); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT); // only create torques and forces for patchy particles

    for (int i = 0; i < runtime; i++)
    {
        cout << i << endl;
        if (i > 0 && i % 20 == 0)
        {
            // cout << "pairs recalculated" << endl;
            delete pairs;
            // delete pairs_onlyb;
            // pairs = obj->calculatepairs(boxes, 3.5);

            // ADD THE PARTICLE ADDITION METHOD HERE
            // Add a random particle to a polymer
            double rand1 = (double)rand()/(double)RAND_MAX;
            // cout << prod << " " << rand1 << endl;
            // pausel();

            if(rand1<prod && indices_to_add.size()>0) {
                //choose a polymer at random
                int polch =  rand() % NpolB;
                int plusorminus = rand() % 2;
                vector1<double> x3(3);
                if(plusorminus == 0) {
                    int index1 = allpolymers[polch][0];
                    int index2 = allpolymers[polch][1];
                    vector1<double> x1 = obj->get_particle(index1);
                    vector1<double> x2 = obj->get_particle(index2);

                    double x11 =  scalar(x1,x1);

                    double x12 = scalar(x1, x2);

                    double x22 = scalar(x2, x2);

                    

                    double t = -(1.1 / sqrt(x11 - 2 * x12 + x22));
                    
                    x3=(1-t)*x1 + t*x2;

                    

                }
                else{
                    int index1 = allpolymers[polch][allpolymers[polch].size()-1];
                    int index2 = allpolymers[polch][allpolymers[polch].size()-2];
                    vector1<double> x1 = obj->get_particle(index1);
                    vector1<double> x2 = obj->get_particle(index2);

                    double x11 = scalar(x1, x1);

                    double x12 = scalar(x1, x2);

                    double x22 = scalar(x2, x2);

                    double t = -(1.1 / sqrt(x11 - 2 * x12 + x22));

                    

                    x3 = (1 - t) * x1 + t * x2;

      
                }

                // vector1<double> randvec(3);
                // randvec[0]=0.01*(-1+(2.*rand())/(double)(RAND_MAX));
                // randvec[1] = 0.01 * (-1 + (2. * rand()) / (double)(RAND_MAX));
                // randvec[2] = 0.01 * (-1 + (2. * rand()) / (double)(RAND_MAX));

                // x3 += randvec;



                    bool overlap = false; 
                    for(int ij = 0  ; ij < indices_everything.size() ; ij++) {
                    double dis = obj->getgeo().distance(obj->get_particle(indices_everything[ij]),x3);
                    if(indices_everything[ij]<totnp) {
                        if (dis < 1.5)
                        {

                            overlap = true;
                            goto ovf;
                        }
                    }
                    else{
                    if(dis<1.09) {
    
                        overlap = true; goto ovf;}
                    }
                    }
                    ovf:
                    if(overlap) {
                        //do not add the particle
                    }
                    else{
                        int ij = rand() % indices_to_add.size();
                        int inde = indices_to_add[ij];

                        // vector1<double> kx(indices_everything.size());
                        // for(int ij = 0  ; ij < indices_everything.size() ; ij++) {
                        //     double dis = obj->getgeo().distance(obj->get_particle(indices_everything[ij]),x3);
                        //     kx[ij] = dis;

                        // }
                        // cout << x3 << endl;
                        // cout << overlap << endl;
                        // cout << kx << endl;
                        // cout << minval(kx) << endl;
                        // cout << minindex(kx) << endl;

                        // pausel();

                        remove_at(indices_to_add, ij);
                        indices_everything.push_back(inde);

                        if(plusorminus == 0) {
                            allpolymers[polch].push_back(inde);
                            std::rotate(allpolymers[polch].rbegin(), allpolymers[polch].rbegin() + 1, allpolymers[polch].rend());
                            mdtriplet tri(allpolymers[polch][0], allpolymers[polch][1], allpolymers[polch][2]);
                            mdpair par(allpolymers[polch][0], allpolymers[polch][1]);
                            polpairs.push_back(par);
                            triplets.push_back(tri);
                        }
                        else{
                            allpolymers[polch].push_back(inde);
                            int nn = allpolymers[polch].size();
                            mdtriplet tri(allpolymers[polch][nn-3], allpolymers[polch][nn-2], allpolymers[polch][nn-1]);
                            mdpair par(allpolymers[polch][nn-2], allpolymers[polch][nn-1]);
                            polpairs.push_back(par);
                            triplets.push_back(tri);
                        }

                        obj->set_particle(x3, inde);
                    }

                    
                }
            


            

            pairs = obj->calculatepairs_parallel(boxes, indices_everything, 3.5);

            // pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_patchy, 3.5);
        }

        obj->advancemom_halfstep(F, T, indices_everything);
        // cout << "advance mom" << endl;
        obj->advance_pos(indices_everything);
        // cout << "advnace pos" << endl;
        // obj->rotate(indices_patchy); // only update the patchy particles with the rotate algorithm

        // cout << "begin forces" << endl;
        F = obj->calculateforces(*pairs, wsa); // calculate the forces due to hard sphere forces

        // cout << "f1" << endl;
        F += obj->calculateforcesharmonic(bindingpairs, diss, myshell.k);

        // cout << "f2" << endl;
        F += obj->calculateforcesdelauny(quads, myshell.kappa);

        // cout << "f3" << endl;
        F += obj->calculateforces(polpairs, pol);

        // cout << "f4" << endl;
        F += obj->calculateforces_threebody(triplets, bp);

        // cout << "f5" << endl;
        
        T.reset(0.0);


        // obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

        generate_uniform_random_matrix(RT); // only generate random torques for the patchy particles
        // cout << "gen ran" << endl;

        obj->create_forces_and_torques_sphere(F, T, RT); // only create torques and forces for patchy particles
        // cout << "create forces" << endl;

        obj->advancemom_halfstep(F, T, indices_everything);
        // cout << "advance mom" << endl;
        if (i % every == 0 && i > 0)
        {

            // cout << i << endl;

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
            // elli += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();
            elli += ss.str();

            poss += extension;
            oris += extension;
            elli += extension;

            string orie = "orient";
            orie += extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            ofstream myfile3;
            myfile3.open(elli.c_str(), std::ios_base::app);

            ofstream myfile4;
            myfile4.open(orie.c_str());

            myfile <<= pos;
            for (int ik = 0; ik < indices_everything.size(); ik++)
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
    // delete pairs_onlyb;
}

#endif