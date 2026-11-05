#ifndef NANOTUBE4_CPP
#define NANOTUBE4_CPP

GeneralPatch CreateGeneralPatchFan(double int1, double size, double range, double ang, geneticcode &g, int Nf)
{

    // in this function the patches are incorporated within a separate structure

    int nt = g.no_types ; // the +1 corresponds to the other material

    vector1<int> vec1(nt);
    for (int i = 0; i < nt; i++)
        vec1[i] = (*(g.patch_num))[i];

    vector1<int> numb(nt);

    for (int i = 0; i < nt; i++)
        numb[i] = Nf + (i + 1) * 5000 + 1;

    int tot = 0;

    for (int i = 0; i < nt; i++)
    {
        for (int j = i; j < nt; j++)
        {
            tot += vec1[i] * vec1[j];
        }
    }

    matrix<double> params(tot, 3);
    int iter = 0;
    int itr = 0;

    for (int i = 0; i < nt; i++)
    {
        for (int j = i; j < nt; j++)
        {
            for (int k = 0; k < g.interactions[itr].size(); k++)
            {
                params(iter, 0) = g.interactions[itr][k] * int1;
                params(iter, 1) = range * size;
                params(iter, 2) = ang;
                iter++;
            }
            itr++;
        }
    }

    int nott = 0;
    for (int i = 0; i < vec1.getsize(); i++)
        nott += vec1[i];

    matrix<double> orient(nott, 3);

    matrix<double> orie(6, 3);

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

    for (int i = 0; i < g.no_types; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            if ((*(g.patch_pos))(i, j) == 1)
            {

                orient(iter2, 0) = orie(j, 0);
                orient(iter2, 1) = orie(j, 1);
                orient(iter2, 2) = orie(j, 2);
                iter2++;
            }
        }
    }

    GeneralPatch c(vec1, numb, params, orient);
    return c;
}



void NanotubeAssembly::run_fan_equil(int runtime, int every, geneticcode &g, ShellProperties &myfan, string strbase)
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

    int Nfan = myfan.posi.getnrows();

    int NN = obj->getN() + Nfan;

    // we also need to define the position of the "particle" that acts as a wall on the

    matrix<double> newdat(NN, 3);



    for(int i = 0 ; i < Nfan ; i++) {
        for(int j = 0 ; j < 3 ; j++) {
        newdat(i,j) = myfan.posi(i,j);
        }
    }
    newdat(Nfan, 0) = 1.;
    newdat(Nfan, 1) = 1.;
    newdat(Nfan, 2) = 1.;

    // newdat(0, 0) = ll / 2.;
    // newdat(0, 1) = ll / 2.;
    // newdat(0, 2) = h;

    obj->initialize(newdat);

    // start with a single particle
    vector<vector<int>> indices(no_types);
    for (int i = 0; i < no_types; i++)
    {
        vector<int> a;
        indices[i] = a;
    }

    indices[0].push_back(Nfan); // add 1 particle to start

    vector<vector<int>> indices_to_add(no_types);
    for (int j = Nfan+1; j < 5000+Nfan; j++)
        indices_to_add[0].push_back(j);

    for (int i = 1; i < no_types; i++)
    {
        vector<int> a;
        indices_to_add[i] = a;
        for (int j = Nfan + i * 5000; j < Nfan+ (i + 1) * 5000; j++)
            indices_to_add[i].push_back(j);
    }

    vector1<int> counts(no_types + 1);
    counts[0] = Nfan;    // index start type 0
    counts[1] = 4999+Nfan; // index start type 1
    for (int j = 2; j < no_types + 1; j++)
    {
        counts[j] = Nfan + counts[j - 1] + 5000;
    }

    vector<int> indices_fan(Nfan);
    for(int i = 0  ; i < Nfan ; i++) {
        indices_fan[i]=i;
    }

    vector<int> indices_combine = flatten(indices);


    vector<int> index_all;
    index_all.reserve(indices_fan.size() + indices_combine.size()); // preallocate memory
    index_all.insert(index_all.end(), indices_fan.begin(), indices_fan.end());
    index_all.insert(index_all.end(), indices_combine.begin(), indices_combine.end());
    ///indices_combine.insert(indices_combine.end(),indices_fan.begin(),indices_fan.end());



    matrix<int> *pairs_onlyhs = obj->calculatepairs_parallel(boxes, index_all, 3.5); // hard sphere

    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_combine, 3.5); //hard sphere
   
    matrix<int> pairs_fan =  myfan.par;




    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);


    WCAPotential wsa(3.0, 1.0, 0.0);
    FENEPotential fene(30.0,1.5);




    F = obj->calculateforcesharmonic(pairs_fan,myfan.k);


    F += obj->calculateforces(*pairs_onlyhs, wsa); // calculate the forces due to hard sphere forces



    for(int j = 0  ; j < myfan.index_anchored.size() ; j++) {
        vector1<double> dFtemp = myfan.anchors[j].force(obj->getcoordinate(myfan.index_anchored[j], 0), obj->getcoordinate(myfan.index_anchored[j], 1), obj->getcoordinate(myfan.index_anchored[j], 2) );
        for (int k = 0; k < 3; k++)
            F(myfan.index_anchored[j], k) += dFtemp[k];
    }



    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

  

    generate_uniform_random_matrix(RT, index_all); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT, index_all, false); // only create torques and forces for patchy particles


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


        obj->advancemom_halfstep(F, T, index_all);


        obj->advance_pos(index_all);

 
        obj->rotate(indices_combine);



        if (i % MC_Move == 0)
        {
            int type_choice = rand() % no_types;
            int ins = rand() % 2;

            if (ins == 0)
            { // insert a particle
                double rangex = ll;
                double rangey = ll;
                double rangez = ll;
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

                for(int j = 0 ; j < indices_fan.size() ; j++)  {
                    double dis = obj->distance(myindex, indices_fan[j]);
                    if (dis < 2.)
                    {
                        en += wsa.energy(dis);
                    }
                    if (en > 100.)
                    {
                        goto energ;
                    }
                }

                for (int j = 0; j < indices_combine.size(); j++)
                {

                    double dis = obj->distance(myindex, indices_combine[j]);
                    if (dis < 2.)
                    {
                        en += obj->particle_energy(myindex, indices_combine[j], *pots);
                        en += wsa.energy(dis);
                    }
                    if (en > 100.)
                    {
                        goto energ;
                    }
                    // if any overlaps are found the particle will be rejected and we can skip the rest
                }
            energ:

                double V = rangex * rangey * rangez;

                double mymu = partialmus[type_choice];
                // cout << "add: " <<  myindex << " " << mymu << " " << en << endl;
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
                    for (int j = 0; j < indices_fan.size(); j++)
                    {
                        if (myindex != indices_fan[j])
                        {
                            double dis = obj->distance(myindex, indices_fan[j]);
                            if (dis < 2.)
                            {
                                en += wsa.energy(dis);
                            }
                        }
                    }
                    for (int j = 0; j < indices_combine.size(); j++)
                    {
                        if (myindex != indices_combine[j])
                        {
                            double dis = obj->distance(myindex, indices_combine[j]);
                            if (dis < 2.)
                            {
                                en += wsa.energy(dis);
                                en += obj->particle_energy(myindex, indices_combine[j], *pots);
                            }
                        }
                        // if any overlaps are found the particle will be rejected and we can skip the rest
                    }
                    // that's the energy currently, we want deltaU, so the energy if it's taken away

                    double V = rangex * rangey * rangez;

                    double mymu = partialmus[type_choice];
                    // cout << "remove: " << myindex << " " << mymu << " " << en << endl;

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
        cout << "MC move" << endl;

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


            delete pairs_onlyb;
            delete pairs_onlyhs;

            index_all.clear();
            index_all.reserve(indices_fan.size() + indices_combine.size()); // preallocate memory
            index_all.insert(index_all.end(), indices_fan.begin(), indices_fan.end());
            index_all.insert(index_all.end(), indices_combine.begin(), indices_combine.end());

            pairs_onlyhs = obj->calculatepairs_parallel(boxes, index_all, 3.5);
            pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_combine, 3.5);
            // cout << "done"  << endl;
        }

        cout << "pairs" << endl;
        matrix<double> bondlengths(pairs_fan.getnrows(),1,1.);

        F = obj->calculateforcesharmonic(pairs_fan, bondlengths, myfan.k);

        F += obj->calculateforces(*pairs_onlyhs, wsa); // calculate the forces due to hard sphere forces

        for (int j = 0; j < myfan.index_anchored.size(); j++)
        {
            vector1<double> dFtemp = myfan.anchors[j].force(obj->getcoordinate(myfan.index_anchored[j], 0), obj->getcoordinate(myfan.index_anchored[j], 1), obj->getcoordinate(myfan.index_anchored[j], 2));
            for (int k = 0; k < 3; k++)
                F(myfan.index_anchored[j], k) += dFtemp[k];
        }

        // pausel();
        T.reset(0.0);
        // cout << *pairs_onlyb << endl;

        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

        cout << "F4" << endl;

        generate_uniform_random_matrix(RT, index_all); // only generate random torques for the patchy particles

        cout << "random" << endl;

        obj->create_forces_and_torques_sphere(F, T, RT, index_all, false); // only create torques and forces for patchy particles

        cout << "create forces" << endl;

        obj->advancemom_halfstep(F, T, index_all);
        cout << "adv mom" << endl;

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
            for (int ik = 0; ik < index_all.size(); ik++)
            {
                myfile <<= pos.getrowvector(index_all[ik]);
                myfile << "\n";
                myfile2 << index_all[ik] << endl;
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
    delete pairs_onlyhs;
    

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