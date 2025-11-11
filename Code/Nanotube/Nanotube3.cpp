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


struct anchoringpoints {
    double epsilon;
    double sigma;
    double anchorpointx,anchorpointy,anchorpointz;

    vector1<double> force(double x, double y, double z) {
        double r2 = SQR(x-anchorpointx)+SQR(y-anchorpointy)+SQR(z-anchorpointz);
    
        vector1<double> force(3);
        if(r2 > 4*sigma*sigma) return force;
        else{
        double prefactor=(12*epsilon)/((exp((3*r2)/SQR(sigma))*SQR(1 + exp(-3*r2)/SQR(sigma)))*SQR(sigma));
        force[0] = -prefactor * (x - anchorpointx);
        force[1] = -prefactor * (y - anchorpointy);
        force[2] = -prefactor * (z - anchorpointz);
        return force;
        }
    }


};

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

void NanotubeAssembly::run_box_equil(int runtime, int every, double mass, geneticcode &g, string strbase)
{

    // generate boxes first
    int no_types = g.no_types;
    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    int NN = obj->getN() + 1; // one particle is the wall
    // we also need to define the position of the "particle" that acts as a wall on the

    matrix<double> newdat(NN, 3);
    double hmin = 10.;
    double h = hmin; // the height of our flap to start with

    newdat(0, 0) = ll / 2.;
    newdat(0, 1) = ll / 2.;
    newdat(0, 2) = h;

    obj->initialize(newdat);

    // start with a single particle
    vector<vector<int>> indices(no_types);
    for (int i = 0; i < no_types; i++)
    {
        vector<int> a;
        indices[i] = a;
    }

    indices[0].push_back(1);

    vector<vector<int>> indices_to_add(no_types);
    for (int j = 2; j < 5000; j++)
        indices_to_add[0].push_back(j);

    for (int i = 1; i < no_types; i++)
    {
        vector<int> a;
        indices_to_add[i] = a;
        for (int j = i * 5000; j < (i + 1) * 5000; j++)
            indices_to_add[i].push_back(j);
    }

    vector1<int> counts(no_types + 1);
    counts[0] = 0;    // index start type 0
    counts[1] = 4999; // index start type 1
    for (int j = 2; j < no_types + 1; j++)
    {
        counts[j] = counts[j - 1] + 5000;
    }

    WCAPotential wsa(3.0, 1.0, 0.0);

    vector<int> temp;
    for (int ty = 0; ty < no_types; ty++)
        for (int i1 = 0; i1 < indices[ty].size(); i1++)
        {
            double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(indices[ty][i1], 2);
            if (dis < 4.)
                temp.push_back(indices[ty][i1]);
        }
    
    vector1<int> pairs_lid(temp.size());
    for (int i1 = 0; i1 < temp.size(); i1++)
    {
        pairs_lid[i1] = temp[i1];
    }

    vector<int> indices_combine = flatten(indices);

    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_combine, 3.5);

    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    double forcep1 = 0; // the force on the lid
    double momp1 = 0;
    // double mass = 10.;

    F = obj->calculateforces(*pairs_onlyb, wsa); // calculate the forces due to hard sphere forces

    for (int i1 = 0; i1 < (pairs_lid).getsize(); i1++)
    {
        double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(pairs_lid[i1], 2);
        double F1 = wsa.force(dis);
        forcep1 += F1;             // force pushes lid up
        F(pairs_lid[i1], 2) -= F1; // pushes particle down
    } // random force
    forcep1 += sqrt(2 * ((*obj).getgamma()) * mass * (*obj).getkT() / (*obj).getdt()) * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808); // temperature
    if (h > hmin)
        forcep1 -= mass * (h - hmin); // constant downwards force if above hmin.

    

    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

    generate_uniform_random_matrix(RT, indices_combine); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT, indices_combine, false); // only create torques and forces for patchy particles

    // particle_adder vv;
    // vv.set_indices(indices_to_add);
    // vv.set_counts(counts);

    // vector1<int> weights;

    // weights = *(g.proportions);

    // we want to choose a
    double mu = g.rate; // can potentially set different chemical potentials

    cout << mu << endl;
    // cout << *(g.proportions) << endl;
    // pausel();
    // box_vol vol2;
    // vol2.ll = ll;
    // vol2.llx = ll;
    // vol2.lly = ll;
    // vol2.llz = obj->getcoordinate(0, 2);
    // vv.set_volume(vol2);
    // vv.set_rate(g.rate);
    // vv.set_irreducible_weights(weights);

    //we choose mu such that under an ideal gas, the total number of particles will be 
    // <N> =  V exp(\mu) where V is the size of our chamber. If we choose individual \mu_i for each
    // particle then we have to choose individual mu. However, we need to solve it such that it still follows
    // the ideal gas scaling for the *total* number of particles, as this is going to be the most relevant
    // quantity in pushing the lid, that we want to keep under control. 

    // so we have exp(p1 x) + exp(p2 x) up to the total number of types and we want it to be equal to exp(mu)


    //perhaps this is not what we want though, we want the reserve to be larger? This would be a bit like playing with the 
    // yes I think what I will do is choose particles non randomly according to the proportion;

    vector1<double> partialmus(no_types);

    if(no_types == 1) {
        //proportions do not matter
        partialmus[0] =  mu;
    }
    else if(no_types == 2) {    
        int xp1 = (*(g.proportions))[0];
        int xp2 = (*(g.proportions))[1];
        // if(xp1==xp2) { //we do not need to do anything
        //     partialmus[0] = mu;
        //     partialmus[1] = mu;
        // }
        // else{
            /*
            this code sets proportions of chemical potential, the bottom one sets proportions of abundances

        double mu2 = 0.0;
        for(int k = 0 ; k < 10 ; k++) {
            //mu2 = mu2 -  ( ( exp(xp1*mu2)+exp(xp2*mu2) -exp(mu) )/(xp1*exp(xp1*mu2)+xp2*(exp(xp2*mu2)) ) );
            mu2 = mu2 - (log((exp(xp1 * mu2) + exp(xp2 * mu2) )) - mu) / ((xp1 * exp(xp1 * mu2) + xp2 * (exp(xp2 * mu2) )) / (exp(xp1 * mu2) + exp(xp2 * mu2) ));
        }
        partialmus[0]=xp1*mu2;
        partialmus[1]=xp2*mu2;
            */
            double mu2 = log(exp(mu) / ((double)(xp1 + xp2 )));
            partialmus[0] = log(double(xp1)) + mu2;
            partialmus[1] = log(double(xp2)) + mu2;
        //}
    }
    else if(no_types == 3) {
        int xp1 = (*(g.proportions))[0];
        int xp2 = (*(g.proportions))[1];
        int xp3 = (*(g.proportions))[2];
        // if(xp1 ==  xp2 && xp2 == xp3) {
        //     partialmus[0] = mu;
        //     partialmus[1] = mu;
        //     partialmus[2] = mu;
        // }
        // else{
            /*
            double mu2 = 0.;
            for (int k = 0; k < 10; k++)
            {
                mu2 = mu2 - (log((exp(xp1 * mu2) + exp(xp2 * mu2) + exp(xp3 * mu2))) - mu) / ((xp1 * exp(xp1 * mu2) + xp2 * (exp(xp2 * mu2) + xp3 * (exp(xp3 * mu2)))) / (exp(xp1 * mu2) + exp(xp2 * mu2) + exp(xp3 * mu2)));
                //cout << mu2 << endl;
            }
            partialmus[0] = xp1 * mu2;
            partialmus[1] = xp2 * mu2;
            partialmus[2] = xp3 * mu2;
            */
           double mu2 = log(exp(mu)/((double)(xp1+xp2+xp3)));
           partialmus[0] = log(double(xp1)) + mu2;
           partialmus[1] = log(double(xp2)) + mu2;
           partialmus[2] = log(double(xp3)) + mu2;
        

    }
    else{

    }


    int MC_Move = 200;

    static vector1<double> total_time(8);


    for (int i = 0; i < runtime; i++)
    {
        //cout << i << endl;
        auto start = std::chrono::high_resolution_clock::now();

        if(i % MC_Move == 0 ) {
            int type_choice = rand() % no_types;
            int ins = rand() % 2;


            if(ins == 0 ) { //insert a particle
                double rangex = ll;
                double rangey = ll;
                double rangez = obj->getcoordinate(0, 2);
                double r1 = rangex * ((double)rand() / (double)(RAND_MAX));
                double r2 = rangey * ((double)rand() / (double)(RAND_MAX));
                double r3 = rangez * ((double)rand() / (double)(RAND_MAX));

                vector1<double> myvec1(3);
                myvec1[0] = r1;
                myvec1[1] = r2;
                myvec1[2] = r3;

                //we can set the particle to be at a particular place first, then calculate the energies
                int index_which = rand() % indices_to_add[type_choice].size();
                int myindex = indices_to_add[type_choice][index_which];

                obj->set_particle(myvec1, myindex);
                double NNN = (double)indices[type_choice].size();

                double en = 0;
                en += wsa.energy(rangez-r3);


                 for(int j = 0 ; j < indices_combine.size() ; j++) {

                     double dis = obj->distance(myindex,indices_combine[j]);
                    if(dis < 2.) { en += wsa.energy(dis);
                    en += obj->particle_energy(myindex, indices_combine[j], *pots);
                    }
                    if(en > 1000. ) {goto energ; }
                    //if any overlaps are found the particle will be rejected and we can skip the rest
                 }
                    energ:

                 double V = rangex*rangey*rangez;

                 double mymu = partialmus[type_choice];
                 //cout << "add: " <<  myindex << " " << mymu << " " << en << endl;
                 double acceptance = min(1.,(V/(1.+NNN))*exp(mymu)*exp(-en));

                 double r = (double) rand()/(double)RAND_MAX;

                 if(r < acceptance) {
                     indices[type_choice].push_back(myindex); //particle is added
                     remove_at(indices_to_add[type_choice], index_which); 
                 }
                 else{

                 }


            }
            else{ //remove a particle
                //choose a random particle;
                if(indices[type_choice].size() > 0 ) {
                    int index_which = (rand() % indices[type_choice].size());

                    int myindex =  indices[type_choice][index_which];

                    double NNN = (double)indices[type_choice].size();

                    double rangex = ll;
                    double rangey = ll;
                    double rangez = obj->getcoordinate(0, 2);

                    double en = 0.0;
                    en += wsa.energy(rangez - obj->getcoordinate(myindex,2));

                    for (int j = 0; j < indices_combine.size() ; j++)
                    {
                        if(myindex != indices_combine[j]) {
                        double dis = obj->distance(myindex, indices_combine[j]);
                        if (dis < 2.) {
                            en += wsa.energy(dis);
                        en += obj->particle_energy(myindex, indices_combine[j], *pots);
                        }
                    }
                        // if any overlaps are found the particle will be rejected and we can skip the rest
                    }
                    //that's the energy currently, we want deltaU, so the energy if it's taken away


                    double V = rangex * rangey * rangez;

                    double mymu = partialmus[type_choice];
                    //cout << "remove: " << myindex << " " << mymu << " " << en << endl;

                    double acceptance = min(1., (NNN / V) * exp(-mymu) * exp(en));

                    double r = (double)rand() / (double)RAND_MAX;

                    if(r < acceptance) {
                        remove_at(indices[type_choice],index_which);

                        auto pos = std::lower_bound(indices_to_add[type_choice].begin(), indices_to_add[type_choice].end(), myindex);
                        indices_to_add[type_choice].insert(pos, myindex);
                    }
                }
            }

            indices_combine = flatten(indices); // all present particles updated

        }

        auto end = std::chrono::high_resolution_clock::now();
        total_time[0] += std::chrono::duration<double>(end - start).count();

        //after the particle is added, pairs need to be recomputed. This happens so long as i % MC and i % pairs are both zero
        start = std::chrono::high_resolution_clock::now();

        if (i > 0 && i % 20 == 0)
        {

            // pairs = obj->calculatepairs(boxes, 3.5);

            // ADD THE PARTICLE ADDITION METHOD HERE
            // bool dd = false;
            // vector1<double> v1(3);
            // int fi;

            // cout << obj->getcoordinate(0, 2) << endl;
            // box_vol vol2;
            // vol2.ll = ll;
            // vol2.llx = ll;
            // vol2.lly = ll;
            // vol2.llz = obj->getcoordinate(0, 2) - 1;

            // // cout << "created box vol" << endl;

            // vv.set_volume(vol2);
            // vv.add_p_w(*obj, indices, dd, v1, fi);
            // cout << "added particle" << endl;
            // if (dd)
            // {

            //     // cout << "added: " << fi << endl;

            //     indices.push_back(fi);

            //     obj->set_particle(v1, fi);
            // }

            vector<int> temp;
            for (int i1 = 0; i1 < indices_combine.size(); i1++)
            {
                double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(indices_combine[i1], 2);
                if (dis < 4.)
                    temp.push_back(indices_combine[i1]);
            }

            // delete pairs_lid;
            pairs_lid.resize(temp.size());
            for (int i1 = 0; i1 < temp.size(); i1++)
            {
                pairs_lid[i1] = temp[i1];
            }
            // cout << "checking lids" << endl;

            delete pairs_onlyb;

            pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_combine, 3.5);
            // cout << "done"  << endl;
        }
        end = std::chrono::high_resolution_clock::now();
        total_time[1] += std::chrono::duration<double>(end - start).count();

        start = std::chrono::high_resolution_clock::now();
        momp1 = (1 - 0.5 * (*obj).getdt() * ((*obj).getgamma()) / mass) * momp1 + ((*obj).getdt() / 2.) * forcep1;
        obj->advancemom_halfstep(F, T, indices_combine);
        obj->advance_pos(indices_combine);
        end = std::chrono::high_resolution_clock::now();
        total_time[2] += std::chrono::duration<double>(end - start).count();

        h = h + ((*obj).getdt() / mass) * momp1; // how heavy do we make the wall? Choose mass = 100.;
        if (h < hmin)
            h = hmin + (hmin - h);   // reflect if hit the bottom
        obj->setcoordinate(0, 2, h); // update in storage

        start = std::chrono::high_resolution_clock::now();
        obj->rotate(indices_combine);
        end = std::chrono::high_resolution_clock::now(); // only update the patchy particles with the rotate algorithm
        total_time[3] += std::chrono::duration<double>(end - start).count();

        start = std::chrono::high_resolution_clock::now();
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

        // stochastic force
        forcep1 += sqrt(2 * ((*obj).getgamma()) * (*obj).getkT() / (*obj).getdt()) * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808); // temperature
        // cout << forcep1 << " ";
        if (h > hmin)
            forcep1 -= mass * (h - hmin); // constant downwards force
        end = std::chrono::high_resolution_clock::now();
        total_time[4] += std::chrono::duration<double>(end - start).count();
        // cout << forcep1 << endl;
        // pausel();
        T.reset(0.0);
        // cout << *pairs_onlyb << endl;
        start = std::chrono::high_resolution_clock::now();
        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy
        end = std::chrono::high_resolution_clock::now();

        total_time[5] += std::chrono::duration<double>(end - start).count();

        start = std::chrono::high_resolution_clock::now();
        generate_uniform_random_matrix(RT,indices_combine); // only generate random torques for the patchy particles
        end = std::chrono::high_resolution_clock::now();
        total_time[6] += std::chrono::duration<double>(end - start).count();

        start = std::chrono::high_resolution_clock::now();
        obj->create_forces_and_torques_sphere(F, T, RT, indices_combine, false); // only create torques and forces for patchy particles
        end = std::chrono::high_resolution_clock::now();
        total_time[7] += std::chrono::duration<double>(end - start).count();


        cout << total_time << endl;

        obj->advancemom_halfstep(F, T, indices_combine);

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
            for (int ik = 0; ik < indices_combine.size(); ik++)
            {
                myfile <<= pos.getrowvector(indices_combine[ik]);
                myfile << "\n";
                myfile2 << indices_combine[ik] << endl;
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

    // we can just add the first element of each
}

void NanotubeAssembly::run_box_equil_cont(int runtime, int every, int startno, double mass, geneticcode &g, string strbase, matrix<double> &pos, matrix<double> &ori, matrix<int> &ind)
{

    // generate boxes first
    int no_types = g.no_types;
    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    int NN = obj->getN() + 1; // one particle is the wall
    // we also need to define the position of the "particle" that acts as a wall on the

    matrix<double> newdat(NN, 3);
    double hmin = 10.;
    double h = hmin; // the height of our flap to start with

    // newdat(0, 0) = ll / 2.;
    // newdat(0, 1) = ll / 2.;
    // newdat(0, 2) = h;


    

    // start with a single particle
    vector<vector<int>> indices(no_types);
    for (int i = 0; i < no_types; i++)
    {
        vector<int> a;
        indices[i] = a;
    }

    //indices[0].push_back(1);
    

    for(int i = 0; i < ind.getnrows() ; i++) {
        int which_type = (int)(ind(i,0)/5000);
        indices[which_type].push_back(ind(i,0));
    }


    vector<vector<int>> indices_to_add(no_types);
    for (int j = 2; j < 5000; j++)
        indices_to_add[0].push_back(j);

    for (int i = 1; i < no_types; i++)
    {
        vector<int> a;
        indices_to_add[i] = a;
        for (int j = i * 5000; j < (i + 1) * 5000; j++)
            indices_to_add[i].push_back(j);
    }


    for(int kk = 0  ; kk < no_types ; kk++) {
        std::unordered_set<int> to_remove(indices[kk].begin(), indices[kk].end());

        indices_to_add[kk].erase(std::remove_if(indices_to_add[kk].begin(), indices_to_add[kk].end(),
                           [&](int x)
                           { return to_remove.count(x); }),
            indices_to_add[kk].end());
    }
    vector<int> indices_combine = flatten(indices);

    newdat(0, 0) = pos(0,0);
    newdat(0, 1) = pos(0,1);
    newdat(0, 2) = pos(0,2);

    h = newdat(0,2);

    for(int j = 1 ; j < pos.getnrows() ; j++) {
        newdat(indices_combine[j-1], 0) = pos(j, 0);
        newdat(indices_combine[j-1], 1) = pos(j, 1);
        newdat(indices_combine[j-1], 2) = pos(j, 2);
    }



    obj->initialize(newdat);



    obj->setorientation(ori);



    WCAPotential wsa(3.0, 1.0, 0.0);

    vector<int> temp;
    for (int ty = 0; ty < no_types; ty++)
        for (int i1 = 0; i1 < indices[ty].size(); i1++)
        {
            double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(indices[ty][i1], 2);
            if (dis < 4.)
                temp.push_back(indices[ty][i1]);
        }

    vector1<int> pairs_lid(temp.size());
    for (int i1 = 0; i1 < temp.size(); i1++)
    {
        pairs_lid[i1] = temp[i1];
    }




    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_combine, 3.5);


    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    double forcep1 = 0; // the force on the lid
    double momp1 = 0;
    // double mass = 10.;

    F = obj->calculateforces(*pairs_onlyb, wsa); // calculate the forces due to hard sphere forces



    for (int i1 = 0; i1 < (pairs_lid).getsize(); i1++)
    {
        double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(pairs_lid[i1], 2);
        double F1 = wsa.force(dis);
        forcep1 += F1;             // force pushes lid up
        F(pairs_lid[i1], 2) -= F1; // pushes particle down
    } // random force
    forcep1 += sqrt(2 * ((*obj).getgamma()) * mass * (*obj).getkT() / (*obj).getdt()) * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808); // temperature
    if (h > hmin)
        forcep1 -= mass * (h - hmin); // constant downwards force if above hmin.



    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy



    generate_uniform_random_matrix(RT, indices_combine); // only generate random torques for the patchy particles



    obj->create_forces_and_torques_sphere(F, T, RT, indices_combine, false); // only create torques and forces for patchy particles



    // particle_adder vv;
    // vv.set_indices(indices_to_add);
    // vv.set_counts(counts);

    // vector1<int> weights;

    // weights = *(g.proportions);

    // we want to choose a
    double mu = g.rate; // can potentially set different chemical potentials

    cout << mu << endl;
    // cout << *(g.proportions) << endl;
    // pausel();
    // box_vol vol2;
    // vol2.ll = ll;
    // vol2.llx = ll;
    // vol2.lly = ll;
    // vol2.llz = obj->getcoordinate(0, 2);
    // vv.set_volume(vol2);
    // vv.set_rate(g.rate);
    // vv.set_irreducible_weights(weights);

    // we choose mu such that under an ideal gas, the total number of particles will be
    //  <N> =  V exp(\mu) where V is the size of our chamber. If we choose individual \mu_i for each
    //  particle then we have to choose individual mu. However, we need to solve it such that it still follows
    //  the ideal gas scaling for the *total* number of particles, as this is going to be the most relevant
    //  quantity in pushing the lid, that we want to keep under control.

    // so we have exp(p1 x) + exp(p2 x) up to the total number of types and we want it to be equal to exp(mu)

    // perhaps this is not what we want though, we want the reserve to be larger? This would be a bit like playing with the
    //  yes I think what I will do is choose particles non randomly according to the proportion;

    vector1<double> partialmus(no_types);

    if (no_types == 1)
    {
        // proportions do not matter
        partialmus[0] = mu;
    }
    else if (no_types == 2)
    {
        int xp1 = (*(g.proportions))[0];
        int xp2 = (*(g.proportions))[1];
        // if(xp1==xp2) { //we do not need to do anything
        //     partialmus[0] = mu;
        //     partialmus[1] = mu;
        // }
        // else{
        /*
        this code sets proportions of chemical potential, the bottom one sets proportions of abundances

    double mu2 = 0.0;
    for(int k = 0 ; k < 10 ; k++) {
        //mu2 = mu2 -  ( ( exp(xp1*mu2)+exp(xp2*mu2) -exp(mu) )/(xp1*exp(xp1*mu2)+xp2*(exp(xp2*mu2)) ) );
        mu2 = mu2 - (log((exp(xp1 * mu2) + exp(xp2 * mu2) )) - mu) / ((xp1 * exp(xp1 * mu2) + xp2 * (exp(xp2 * mu2) )) / (exp(xp1 * mu2) + exp(xp2 * mu2) ));
    }
    partialmus[0]=xp1*mu2;
    partialmus[1]=xp2*mu2;
        */
        double mu2 = log(exp(mu) / ((double)(xp1 + xp2)));
        partialmus[0] = log(double(xp1)) + mu2;
        partialmus[1] = log(double(xp2)) + mu2;
        //}
    }
    else if (no_types == 3)
    {
        int xp1 = (*(g.proportions))[0];
        int xp2 = (*(g.proportions))[1];
        int xp3 = (*(g.proportions))[2];
        // if(xp1 ==  xp2 && xp2 == xp3) {
        //     partialmus[0] = mu;
        //     partialmus[1] = mu;
        //     partialmus[2] = mu;
        // }
        // else{
        /*
        double mu2 = 0.;
        for (int k = 0; k < 10; k++)
        {
            mu2 = mu2 - (log((exp(xp1 * mu2) + exp(xp2 * mu2) + exp(xp3 * mu2))) - mu) / ((xp1 * exp(xp1 * mu2) + xp2 * (exp(xp2 * mu2) + xp3 * (exp(xp3 * mu2)))) / (exp(xp1 * mu2) + exp(xp2 * mu2) + exp(xp3 * mu2)));
            //cout << mu2 << endl;
        }
        partialmus[0] = xp1 * mu2;
        partialmus[1] = xp2 * mu2;
        partialmus[2] = xp3 * mu2;
        */
        double mu2 = log(exp(mu) / ((double)(xp1 + xp2 + xp3)));
        partialmus[0] = log(double(xp1)) + mu2;
        partialmus[1] = log(double(xp2)) + mu2;
        partialmus[2] = log(double(xp3)) + mu2;
    }
    else
    {
    }



    int MC_Move = 200;


    for (int i = 0; i < runtime; i++)
    {
        cout << "time: " << i << endl;
        if (i % MC_Move == 0 )
        {
            int type_choice = rand() % no_types;
            int ins = rand() % 2;

            if (ins == 0)
            { // insert a particle
                double rangex = ll;
                double rangey = ll;
                double rangez = obj->getcoordinate(0, 2);
                double r1 = rangex * ((double)rand() / (double)(RAND_MAX));
                double r2 = rangey * ((double)rand() / (double)(RAND_MAX));
                double r3 = rangez * ((double)rand() / (double)(RAND_MAX));

                vector1<double> myvec1(3);
                myvec1[0] = r1;
                myvec1[1] = r2;
                myvec1[2] = r3;

                // we can set the particle to be at a particular place first, then calculate the energies
                int index_which = rand() % indices_to_add[type_choice].size();
                int myindex = indices_to_add[type_choice][index_which];

                obj->set_particle(myvec1, myindex);
                double NNN = (double)indices[type_choice].size();

                double en = 0;
                en += wsa.energy(rangez - r3);

                for (int j = 0; j < indices_combine.size(); j++)
                {

                    double dis = obj->distance(myindex, indices_combine[j]);
                    if (dis < 2.)
                        en += wsa.energy(dis);
                    en += obj->particle_energy(myindex, indices_combine[j], *pots);
                    if (en > 1000.)
                    {
                        goto energ;
                    }
                    // if any overlaps are found the particle will be rejected and we can skip the rest
                }
            energ:

                double V = rangex * rangey * rangez;

                double mymu = partialmus[type_choice];
                cout << "add: " << myindex << " " << mymu << " " << en << endl;
                double acceptance = min(1., (V / (1. + NNN)) * exp(mymu) * exp(-en));

                double r = (double)rand() / (double)RAND_MAX;

                if (r < acceptance)
                {
                    indices[type_choice].push_back(myindex); // particle is added
                    remove_at(indices_to_add[type_choice], index_which);
                }
                else
                {
                }
            }
            else
            { // remove a particle
                // choose a random particle;
                if (indices[type_choice].size() > 0)
                {
                    int index_which = (rand() % indices[type_choice].size());

                    int myindex = indices[type_choice][index_which];

                    double NNN = (double)indices[type_choice].size();

                    double rangex = ll;
                    double rangey = ll;
                    double rangez = obj->getcoordinate(0, 2);

                    double en = 0.0;
                    en += wsa.energy(rangez - obj->getcoordinate(myindex, 2));

                    for (int j = 0; j < indices_combine.size(); j++)
                    {
                        if (myindex != indices_combine[j])
                        {
                            double dis = obj->distance(myindex, indices_combine[j]);
                            if (dis < 2.)
                                en += wsa.energy(dis);
                            en += obj->particle_energy(myindex, indices_combine[j], *pots);
                        }
                        // if any overlaps are found the particle will be rejected and we can skip the rest
                    }
                    // that's the energy currently, we want deltaU, so the energy if it's taken away

                    double V = rangex * rangey * rangez;

                    double mymu = partialmus[type_choice];
                    cout << "remove: " << myindex << " " << mymu << " " << en << endl;

                    double acceptance = min(1., (NNN / V) * exp(-mymu) * exp(en));

                    double r = (double)rand() / (double)RAND_MAX;

                    if (r < acceptance)
                    {
                        remove_at(indices[type_choice], index_which);

                        auto pos = std::lower_bound(indices_to_add[type_choice].begin(), indices_to_add[type_choice].end(), myindex);
                        indices_to_add[type_choice].insert(pos, myindex);
                    }
                }
            }

            indices_combine = flatten(indices); // all present particles updated
        }

        // after the particle is added, pairs need to be recomputed. This happens so long as i % MC and i % pairs are both zero

        if (i > 0 && i % 20 == 0)
        {

            // pairs = obj->calculatepairs(boxes, 3.5);

            // ADD THE PARTICLE ADDITION METHOD HERE
            // bool dd = false;
            // vector1<double> v1(3);
            // int fi;

            // cout << obj->getcoordinate(0, 2) << endl;
            // box_vol vol2;
            // vol2.ll = ll;
            // vol2.llx = ll;
            // vol2.lly = ll;
            // vol2.llz = obj->getcoordinate(0, 2) - 1;

            // // cout << "created box vol" << endl;

            // vv.set_volume(vol2);
            // vv.add_p_w(*obj, indices, dd, v1, fi);
            // cout << "added particle" << endl;
            // if (dd)
            // {

            //     // cout << "added: " << fi << endl;

            //     indices.push_back(fi);

            //     obj->set_particle(v1, fi);
            // }

            vector<int> temp;
            for (int i1 = 0; i1 < indices_combine.size(); i1++)
            {
                double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(indices_combine[i1], 2);
                if (dis < 4.)
                    temp.push_back(indices_combine[i1]);
            }

            // delete pairs_lid;
            pairs_lid.resize(temp.size());
            for (int i1 = 0; i1 < temp.size(); i1++)
            {
                pairs_lid[i1] = temp[i1];
            }
            // cout << "checking lids" << endl;

            delete pairs_onlyb;

            pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_combine, 3.5);
            // cout << "done"  << endl;
        }

        momp1 = (1 - 0.5 * (*obj).getdt() * ((*obj).getgamma()) / mass) * momp1 + ((*obj).getdt() / 2.) * forcep1;
        obj->advancemom_halfstep(F, T, indices_combine);
        obj->advance_pos(indices_combine);

        h = h + ((*obj).getdt() / mass) * momp1; // how heavy do we make the wall? Choose mass = 100.;
        if (h < hmin)
            h = hmin + (hmin - h);   // reflect if hit the bottom
        obj->setcoordinate(0, 2, h); // update in storage

        obj->rotate(indices_combine);                // only update the patchy particles with the rotate algorithm
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

        // stochastic force
        forcep1 += sqrt(2 * ((*obj).getgamma()) * (*obj).getkT() / (*obj).getdt()) * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808); // temperature
        // cout << forcep1 << " ";
        if (h > hmin)
            forcep1 -= mass * (h - hmin); // constant downwards force

        // cout << forcep1 << endl;
        // pausel();
        T.reset(0.0);
        // cout << *pairs_onlyb << endl;
        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy
        // cout << "huh" << endl;
        


        generate_uniform_random_matrix(RT, indices_combine); // only generate random torques for the patchy particles

        obj->create_forces_and_torques_sphere(F, T, RT, indices_combine,false); // only create torques and forces for patchy particles


        obj->advancemom_halfstep(F, T, indices_combine);

        momp1 = (1 - 0.5 * (*obj).getdt() * ((*obj).getgamma()) / mass) * momp1 + ((*obj).getdt() / 2.) * forcep1;


        if (i % every == 0 && i > 0)
        {

            // cout << i << endl;

            stringstream ss;

            ss << setw(number_of_digits) << setfill('0') << startno + (i / every);

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
            for (int ik = 0; ik < indices_combine.size(); ik++)
            {
                myfile <<= pos.getrowvector(indices_combine[ik]);
                myfile << "\n";
                myfile2 << indices_combine[ik] << endl;
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

    // we can just add the first element of each

    // cout << total_timeMC << endl;
    // cout << total_timePairs << endl;
    // cout << total_timeForce << endl;
}

void NanotubeAssembly::run_box_anchors(int runtime, int every, double mass, int NumA, geneticcode &g, string strbase)
{

    // generate boxes first
    int no_types = g.no_types;
    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    int NN = obj->getN() + 1; // one particle is the wall
    // we also need to define the position of the "particle" that acts as a wall on the

    matrix<double> newdat(NN, 3);
    double hmin = 10.;
    double h = hmin; // the height of our flap to start with

    newdat(0, 0) = ll / 2.;
    newdat(0, 1) = ll / 2.;
    newdat(0, 2) = h;

    obj->initialize(newdat);

    // start with a single particle
    vector<int> indices;
    indices.reserve(NN); // these are all the particles we will be adding

    indices.push_back(1);
    
    double dl = ll/(double)NumA;

    vector<vector1<double> > anchorpoints;
    for(int ix = 0 ; ix < NumA ; ix++) {
        for(int iy = 0 ; iy < NumA ; iy++) {
            
            vector1<double> pos(3);
            pos[0] = (ix+0.5)*dl;
            pos[1] = (iy+0.5)*dl;
            pos[2] = 0.5;
            anchorpoints.push_back(pos);
        }
    }
    // the density is the  number of points


    // vector<double> indices_weights;
    //  vector<int> indices_to_add;
    //  indices_to_add.reserve(NN-2);
    //  for(int i =2 ; i < NN ; i++) {
    //      indices_to_add.push_back(i); // all the ones to add
    //  }

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
    for (int i = 2; i < no_types * 5000; i++)
        indices_to_add.push_back(i);
    // 1 is already present

    vector1<int> counts(no_types + 1);
    counts[0] = 0;    // index start type 0
    counts[1] = 4999; // index start type 1
    for (int j = 2; j < no_types + 1; j++)
    {
        counts[j] = counts[j - 1] + 5000;
    }

    WCAPotential wsa(3.0, 1.0, 0.0);

    vector<int> temp;
    for (int i1 = 0; i1 < indices.size(); i1++)
    {
        double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(indices[i1], 2);
        if (dis < 4.)
            temp.push_back(indices[i1]);
    }
    vector1<int> pairs_lid(temp.size());
    for (int i1 = 0; i1 < temp.size(); i1++)
    {
        pairs_lid[i1] = temp[i1];
    }

    vector<int> temp2;
    for (int i1 = 0; i1 < indices.size(); i1++)
    {
        double dis = obj->getcoordinate(indices[i1], 2);
        if (dis < 4.)
            temp2.push_back(indices[i1]);
    }
    vector1<int> pairs_anchor(temp2.size());
    for (int i1 = 0; i1 < temp2.size(); i1++)
    {
        pairs_anchor[i1] = temp2[i1];
    }

    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, indices, 3.5);

    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);

    double forcep1 = 0; // the force on the lid
    double momp1 = 0;
    // double mass = 10.;

    F = obj->calculateforces(*pairs_onlyb, wsa); // calculate the forces due to hard sphere forces

    for (int i1 = 0; i1 < (pairs_lid).getsize(); i1++)
    {
        double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(pairs_lid[i1], 2);
        double F1 = wsa.force(dis);
        forcep1 += F1;             // force pushes lid up
        F(pairs_lid[i1], 2) -= F1; // pushes particle down
    } // random force
    forcep1 += sqrt(2 * ((*obj).getgamma()) * mass * (*obj).getkT() / (*obj).getdt()) * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808); // temperature
    if (h > hmin)
        forcep1 -= mass * (h - hmin); // constant downwards force if above hmin.

    anchoringpoints apot;
    apot.epsilon = 10.;
    apot.sigma = 1.;

    for (int i1 = 0; i1 < pairs_anchor.getsize(); i1++)
    {
        for (int k = 0; k < anchorpoints.size(); k++)
        {
            double dis = obj->distance(pairs_anchor[i1], anchorpoints[k]);
            if (dis < 2.)
            {
                apot.anchorpointx = anchorpoints[k][0];
                apot.anchorpointy = anchorpoints[k][1];
                apot.anchorpointz = anchorpoints[k][2];
                vector1<double> df = apot.force(obj->getcoordinate(pairs_anchor[i1], 0), obj->getcoordinate(pairs_anchor[i1], 1), obj->getcoordinate(pairs_anchor[i1], 2));
                F(pairs_anchor[i1], 0) += df[0];
                F(pairs_anchor[i1], 1) += df[1];
                F(pairs_anchor[i1], 2) += df[2];

                break;
            }
        }
    }

    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

    generate_uniform_random_matrix(RT, indices); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT, indices, false); // only create torques and forces for patchy particles

    particle_adder vv;
    vv.set_indices(indices_to_add);
    vv.set_counts(counts);

    vector1<int> weights;

    weights = *(g.proportions);

    // we want to choose a

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

            cout << obj->getcoordinate(0, 2) << endl;
            box_vol vol2;
            vol2.ll = ll;
            vol2.llx = ll;
            vol2.lly = ll;
            vol2.llz = obj->getcoordinate(0, 2) - 1;

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
                double dis = obj->getcoordinate(0, 2) - obj->getcoordinate(indices[i1], 2);
                if (dis < 4.)
                    temp.push_back(indices[i1]);
            }

            // delete pairs_lid;
            pairs_lid.resize(temp.size());
            for (int i1 = 0; i1 < temp.size(); i1++)
            {
                pairs_lid[i1] = temp[i1];
            }

            vector<int> temp2;
            for (int i1 = 0; i1 < indices.size(); i1++)
            {
                double dis =obj->getcoordinate(indices[i1], 2);
                if (dis < 4.)
                    temp2.push_back(indices[i1]);
            }

            // delete pairs_lid;
            pairs_anchor.resize(temp2.size());
            for (int i1 = 0; i1 < temp2.size(); i1++)
            {
                pairs_anchor[i1] = temp2[i1];
            }
            // cout << "checking lids" << endl;

            delete pairs_onlyb;

            pairs_onlyb = obj->calculatepairs_parallel(boxes, indices, 3.5);
            // cout << "done"  << endl;
        }

        momp1 = (1 - 0.5 * (*obj).getdt() * ((*obj).getgamma()) / mass) * momp1 + ((*obj).getdt() / 2.) * forcep1;
        obj->advancemom_halfstep(F, T, indices);
        obj->advance_pos(indices);

        h = h + ((*obj).getdt() / mass) * momp1; // how heavy do we make the wall? Choose mass = 100.;
        if (h < hmin)
            h = hmin + (hmin - h);   // reflect if hit the bottom
        obj->setcoordinate(0, 2, h); // update in storage

        obj->rotate(indices);                        // only update the patchy particles with the rotate algorithm
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

        // stochastic force
        forcep1 += sqrt(2 * ((*obj).getgamma()) * (*obj).getkT() / (*obj).getdt()) * (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808); // temperature
        // cout << forcep1 << " ";
        if (h > hmin)
            forcep1 -= mass * (h - hmin); // constant downwards force

        for (int i1 = 0; i1 < pairs_anchor.getsize(); i1++)
        {
            for (int k = 0; k < anchorpoints.size(); k++)
            {
                double dis = obj->distance(pairs_anchor[i1], anchorpoints[k]);
                if (dis < 2.)
                {
                    apot.anchorpointx = anchorpoints[k][0];
                    apot.anchorpointy = anchorpoints[k][1];
                    apot.anchorpointz = anchorpoints[k][2];
                    vector1<double> df = apot.force(obj->getcoordinate(pairs_anchor[i1], 0), obj->getcoordinate(pairs_anchor[i1], 1), obj->getcoordinate(pairs_anchor[i1], 2));
                    F(pairs_anchor[i1], 0) += df[0];
                    F(pairs_anchor[i1], 1) += df[1];
                    F(pairs_anchor[i1], 2) += df[2];

                    //break;
                }
            }
        }

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
            for (int ik = 0; ik < indices.size(); ik++)
            {
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

    // we can just add the first element of each
}

#endif