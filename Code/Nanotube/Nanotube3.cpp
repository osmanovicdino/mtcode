#ifndef NANOTUBE3_CPP
#define NANOTUBE3_CPP


GeneralPatch CreateGeneralPatch(double int1, double size, double range, double ang, geneticcode &g) {
    
    int nt= g.no_types;
    
    vector1<int> vec1(nt);
    for(int i = 0 ; i < nt ; i++)
        vec1[i] = (*(g.patch_num))[i];
   
    vector1<int> numb(nt);

    for(int i = 0 ; i < nt ; i++)
        numb[i] = (i+1)*5000+1;

    int tot = 0;

    for(int i = 0 ; i < nt ; i++) {
        for(int j = i ; j < nt; j++) {
            tot += vec1[i] * vec1[j];
        }
    }
    
    matrix<double> params(tot,3);
    int iter = 0;
    int itr = 0;

    for(int i = 0  ; i < nt; i++) {
        for(int j = i  ; j < nt ; j++) {
            for(int k = 0  ; k < g.interactions[itr].size() ; k++) {
                params(iter,0) = g.interactions[itr][k]*int1;
                params(iter,1) = range*size;
                params(iter,2) = ang;
                iter++;
            }
            itr++;
        }
    }

    int nott = 0;
    for(int i = 0  ; i < vec1.getsize() ; i++)
        nott+=vec1[i];
    


    matrix<double> orient(nott,3);

    matrix<double> orie(6,3);

    orie(0, 0) = 1.0;
    orie(0, 1) = 0.;
    orie(0, 2) = 0.;

    orie(1, 0) = -1.0;
    orie(1, 1) = 0.;
    orie(1, 2) = 0.;

    orie(2, 0) = 0.;
    orie(2, 1) = 1.;
    orie(2, 2) = 0.;

    orie(3, 0) = 0.;
    orie(3, 1) = -1.;
    orie(3, 2) = 0.;

    orie(4, 0) = 0.0;
    orie(4, 1) = 0.;
    orie(4, 2) = 1.;

    orie(5, 0) = 0.;
    orie(5, 1) = 0.;
    orie(5, 2) = -1.;


    int iter2 = 0;

    
    for(int i = 0  ; i < g.no_types ; i++) {
        for(int j = 0 ; j < 6 ; j++) {
            if( (*(g.patch_pos))(i, j) == 1) {

            orient(iter2, 0) = orie(j, 0);
            orient(iter2, 1) = orie(j, 1);
            orient(iter2, 2) = orie(j, 2);
            iter2++;
            }
        }
    }


    GeneralPatch c(vec1,numb,params,orient);
    return c;
}

void NanotubeAssembly::run_box(int runtime, int every, double mass, geneticcode &g, string strbase)
{

    //generate boxes first
    int no_types =  g.no_types;
    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    int NN = obj->getN() + 1; //one particle is the wall
    //we also need to define the position of the "particle" that acts as a wall on the 

    matrix<double> newdat(NN,3);
    double hmin = 10.;
    double h = hmin; //the height of our flap to start with

    newdat(0,0) = ll/2.;
    newdat(0,1) = ll/2.;
    newdat(0,2) = h;

    obj->initialize(newdat);

    //start with a single particle
    vector<int> indices;
    indices.reserve(NN); //these are all the particles we will be adding

    
    indices.push_back(1);

    //vector<double> indices_weights;
    // vector<int> indices_to_add;
    // indices_to_add.reserve(NN-2);
    // for(int i =2 ; i < NN ; i++) {
    //     indices_to_add.push_back(i); // all the ones to add
    // }

    // vector<vector<int> > indices_to_add;

    // int i1 = 0;
    // vector<int> v;
    // for (int k = 1; k < 5000; k++)
    // {
    //     v.push_back(i1);
    //     i1++;
    // }
    // indices_to_add.push_back(v);

    // for(int j = 1 ; j < no_types ; j++) {
    //     vector<int> v;
    //     int i = 0 ;
    //     for(int k = 0  ; k < 5000 ; k++) {
    //         v.push_back(i);
    //         i++;
    //     }
    //     indices_to_add.push_back(v);
    // }
    vector<int> indices_to_add;
    for(int i = 2 ; i < no_types*5000 ; i++)
        indices_to_add.push_back(i);
    // 1 is already present


    vector1<int> counts(no_types+1);
    counts[0] = 0; //index start type 0
    counts[1] = 4999; //index start type 1
    for(int j = 2 ; j < no_types+1 ; j++) {
        counts[j] = counts[j-1] + 5000;
    }



    WCAPotential wsa(3.0, 1.0, 0.0);

   
    vector<int> temp;
    for(int i1 = 0  ; i1 < indices.size() ; i1++) {
        double dis = obj->getcoordinate(0,2) - obj->getcoordinate(indices[i1],2);
        if(dis < 4.) temp.push_back(indices[i1]);
        }
    vector1<int> pairs_lid(temp.size());
    for(int i1  = 0 ; i1 < temp.size() ; i1++) {
        pairs_lid[i1] = temp[i1];
    }
    


    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, indices, 3.5);

    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    

    double forcep1 = 0; //the force on the lid
    double momp1 = 0;
    // double mass = 10.;

    F = obj->calculateforces(*pairs_onlyb, wsa); // calculate the forces due to hard sphere forces

    for(int i1 = 0  ; i1 < (pairs_lid).getsize() ; i1++) {
        double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(pairs_lid[i1], 2);
        double F1 = wsa.force(dis);
        forcep1 += F1; //force pushes lid up
        F(pairs_lid[i1],2) -= F1; //pushes particle down
    } //random force
    forcep1 += sqrt(2 * ((*obj).getgamma() ) *mass * (*obj).getkT() / (*obj).getdt()) * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808); // temperature
    if(h>hmin)
    forcep1 -= mass*(h-hmin); //constant downwards force if above hmin.
    
    
    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

    generate_uniform_random_matrix(RT, indices); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT, indices, false); // only create torques and forces for patchy particles

    particle_adder vv;
    vv.set_indices(indices_to_add);
    vv.set_counts(counts);

    vector1<int> weights;

    weights = *(g.proportions);

    //we want to choose a 

    box_vol vol2;
    vol2.ll = ll;
    vol2.llx = ll;
    vol2.lly = ll;
    vol2.llz = obj->getcoordinate(0, 2);
    vv.set_volume(vol2);
    vv.set_rate(g.rate);
    vv.set_irreducible_weights(weights);




    for (int i = 0; i < runtime; i++)
    {
        cout << i << endl;
        if (i > 0 && i % 20 == 0)
        {


            // pairs = obj->calculatepairs(boxes, 3.5);
 
            // ADD THE PARTICLE ADDITION METHOD HERE
            bool dd = false;
            vector1<double> v1(3);
            int fi;

            cout << obj->getcoordinate(0,2) << endl;
           box_vol vol2;
            vol2.ll = ll;
            vol2.llx = ll;
            vol2.lly = ll;
            vol2.llz = obj->getcoordinate(0,2)-1;

            // cout << "created box vol" << endl;

            vv.set_volume(vol2);
            vv.add_p_w(*obj, indices, dd, v1, fi);
            // cout << "added particle" << endl;
            if (dd)
            {

                // cout << "added: " << fi << endl;

                indices.push_back(fi);

                obj->set_particle(v1, fi);
            }


            vector<int> temp;
            for (int i1 = 0; i1 < indices.size(); i1++)
            {
                double dis = obj->getcoordinate(0,2) - obj->getcoordinate(indices[i1],2);
                if (dis < 4.)
                    temp.push_back(indices[i1]);
            }

            // delete pairs_lid;
            pairs_lid.resize(temp.size());
            for (int i1 = 0; i1 < temp.size(); i1++)
            {
                pairs_lid[i1] = temp[i1];
            }
            // cout << "checking lids" << endl;

            delete pairs_onlyb;

            pairs_onlyb = obj->calculatepairs_parallel(boxes, indices, 3.5);
            // cout << "done"  << endl;
        }

        momp1 = (1 - 0.5 * (*obj).getdt()*((*obj).getgamma()) / mass) * momp1 + ((*obj).getdt() / 2.) * forcep1;
        obj->advancemom_halfstep(F, T, indices);
        obj->advance_pos(indices);
        
        h = h + ((*obj).getdt() / mass) * momp1; // how heavy do we make the wall? Choose mass = 100.;
        if( h < hmin) h = hmin+(hmin-h); //reflect if hit the bottom
        obj->setcoordinate(0,2,h); //update in storage
        
        obj->rotate(indices); // only update the patchy particles with the rotate algorithm
        F = obj->calculateforces(*pairs_onlyb, wsa); // calculate the forces due to hard sphere forces
        forcep1 = 0;
        for (int i1 = 0; i1 < (pairs_lid).getsize(); i1++)
        {

            double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(pairs_lid[i1], 2);

            double F1 = wsa.force(dis);
            forcep1 += F1;                 // force pushes lid up
            F(((pairs_lid))[i1], 2) -= F1; // pushes particle down
        }
        // cout << forcep1 << " ";

        //stochastic force
        forcep1 += sqrt(2 * ((*obj).getgamma())  * (*obj).getkT() / (*obj).getdt()) * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808); // temperature
        // cout << forcep1 << " ";
        if (h > hmin)
            forcep1 -= mass*(h-hmin); // constant downwards force

        // cout << forcep1 << endl;
        // pausel();
        T.reset(0.0);
        // cout << *pairs_onlyb << endl;
        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy
        // cout << "huh" << endl;
        generate_uniform_random_matrix(RT); // only generate random torques for the patchy particles

        obj->create_forces_and_torques_sphere(F, T, RT); // only create torques and forces for patchy particles
        obj->advancemom_halfstep(F, T, indices);

        momp1 = (1 - 0.5 * (*obj).getdt() * ((*obj).getgamma()) / mass) * momp1 + ((*obj).getdt() / 2.) * forcep1;

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

            myfile <<= pos.getrowvector(0);
            myfile << "\n";
            for (int ik = 0; ik < indices.size(); ik++) {
                myfile <<= pos.getrowvector(indices[ik]);
                myfile << "\n";
                myfile2 << indices[ik] << endl;
            }


            myfile4 <<= obj->getorientation();

            myfile.close();
            myfile2.close();
            myfile3.close();
            myfile4.close();

            // pausel();
        }
    }
    delete pairs_onlyb;

    // in this example, we probably want to specify the weights to be such that the total rate of each remains fixed?

    // i.e. if we have  rates as (rate*p1) and (rate*p2) and rate*p3 where p1 are the proportions we want to keep the ps proportional

    // but the complication is that each weight has to change

    // we can make all the weights point to the same reference? And then as the particle is released we update the reference?

    // we can also make it so that the length of indicies to add is the number of types * 5000 (for instance)

    // there's probably a better way of doing this

    // differential rates

    //we can just add the first element of each
}

#endif