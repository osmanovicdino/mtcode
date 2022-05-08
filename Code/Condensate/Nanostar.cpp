#ifndef NANOSTAR_CPP
#define NANOSTAR_CPP

Nanostar::Nanostar(int N, double ll, int length_of_branchh) :  bindpairs(vector<mdpair>()), bendpairs(vector<mdtriplet>()), stickerList(vector<int>()){
    obj = new LangevinNVT;
    dimension = 3;
    l = ll;

    num_nanostars = N;


    total_particles = 0;

    num_branches = 3;

    length_of_branch = length_of_branchh;

    total_parts_per_nanostar = num_branches*length_of_branch+1;

    total_particles = 0;

    for (int i = 0; i < num_nanostars + 1; i++)
    {
        for(int j = 0  ; j < num_branches ; j++) {
        stickerList.push_back(i*total_parts_per_nanostar + length_of_branch * j);
        }
    }

    vector1<bool> is_periodic(dimension,true);
    cube bc(ll,is_periodic,dimension);

    double sigma = 1.0;

    double att_epp = 10.0;

    WCAPotential StickerPotential(10.0,sigma,att_epp);

    WCAPotential HS(20.0, sigma, 0.0);

    FENEPotential nrr(50,1.4);

    double preferred_angle = 0.0;

    double bending_strength = 10.0;

    BendingPotential nfr2(bending_strength,preferred_angle);



    for(int i =  0 ; i < num_nanostars; i++) {
        this->create_nanostar();
    }

    LangevinNVT b(bc);

    double kT = 1.0;
    double dt = 0.005;
    double eta = 5.;
    // gamma = eta;

    
    // b.setinteractions(nswca);
    b.setkT(kT);
    b.setdt(dt);
    b.setgamma(eta);
    b.setm(1.0);


    *obj = b;

    matrix<double> store = create_initial_state();


    potential *q1 = StickerPotential.clone();
    faa = q1;
    potential *q2 = nrr.clone();
    bindp = q2;
    potential3 *q3 = nfr2.clone();
    bendp = q3;
    potential *q4 = HS.clone();
    hs = q4;

    matrix<double> moms(store.getnrows(), dimension);
    obj->setdat(store);
    obj->setmom(moms);

    

}

matrix<double> rotationmatrix(double theta, double phi, double psi) {
matrix<double> rot(3,3);

rot(0, 0) = cos(theta)*cos(psi);
rot(0, 1) = -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
rot(0, 2) = sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
rot(1, 0) = cos(theta)*sin(psi);
rot(1, 1) = cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);
rot(1, 2) = -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
rot(2, 0) = -sin(theta);
rot(2, 1) = sin(phi)*cos(theta);
rot(2, 2) = cos(phi)*cos(theta);

return rot;

}

matrix<double> Nanostar::create_initial_state() {

    matrix<double> store(total_particles,3);

    vector1<double> ori1(3);
    vector1<double> ori2(3);
    vector1<double> ori3(3);

    double orig1[3] = {0.942809, 0, 0.};
    double orig2[3] = {-0.471405, 0.816497, 0.};
    double orig3[3] = {-0.471405, -0.816497, 0.};

    for (int j = 0; j < 3; j++)
    {
        ori1[j] = 1.1 * orig1[j];
        ori2[j] = 1.1 * orig2[j];
        ori3[j] = 1.1 * orig3[j];
    }



    for(int i = 0  ; i < num_nanostars ; i++) {
        //cout << i << endl;
        double x = l * (double)rand()/(double)(RAND_MAX);
        double y = l * (double)rand() / (double)(RAND_MAX);
        double z = l * (double)rand() / (double)(RAND_MAX);

        vector1<double> centre(3);
        centre[0] = x;
        centre[1] = y;
        centre[2] = z;

        store(i * total_parts_per_nanostar, 0) = x;
        store(i * total_parts_per_nanostar, 1) = y;
        store(i * total_parts_per_nanostar, 2) = z;

        double theta = 2 * pi * ((double)rand() / ((double)RAND_MAX));
        double phi = 2 * pi * ((double)rand() / ((double)RAND_MAX));
        double psi = 2 * pi * ((double)rand() / ((double)RAND_MAX));

        matrix<double> rot(rotationmatrix(theta,phi,psi));
        //do a random rotation on the nanostar
        vector1<double> nori1 = rot * ori1;
        vector1<double> nori2 = rot * ori2;
        vector1<double> nori3 = rot * ori3;

        for (int j = 0; j < length_of_branch; j++)
        {
            vector1<double> pos1 = centre + (double)(j+1)*nori1;
            vector1<double> pos2 = centre + (double)(j+1)*nori2;
            vector1<double> pos3 = centre + (double)(j+1)*nori3;


            obj->getgeo().correct_position(pos1);
            obj->getgeo().correct_position(pos2);
            obj->getgeo().correct_position(pos3);

            for(int k = 0  ; k < dimension ; k++) {
            store(i * total_parts_per_nanostar + 1+0 * (length_of_branch) + j, k) = pos1[k];
            store(i * total_parts_per_nanostar + 1+1 * (length_of_branch) + j, k) = pos2[k];
            store(i * total_parts_per_nanostar + 1+2 * (length_of_branch) + j, k) = pos3[k];
            }
        }
    }
    cout << "nanostars added" << endl;

    return store;

}


matrix<double> Nanostar::create_initial_state(string s)
{

    double T;
    bool err;
    matrix<double> store = importcsv(s,T,err);
    if (store.getNsafe() !=  total_particles)
        error("size of data in filename A not correct (num of particles)");

    if (store.getncols() != dimension)
        error("size of data in filename A not correct (dimension of space)");

    if(err) error("IMPORT FILE NOT FOUND!");

    return store;
}


void Nanostar::set_initial_state(string s)
{

    double T;
    bool err;
    matrix<double> store = importcsv(s, T, err);
    if (store.getNsafe() != total_particles)
        error("size of data in filename A not correct (num of particles)");

    if (store.getncols() != dimension)
        error("size of data in filename A not correct (dimension of space)");

    if (err)
        error("IMPORT FILE NOT FOUND!");

    (*obj).setdat(store);
}

void Nanostar::create_nanostar() {


    int initial = total_particles;

    for(int i = 0 ; i < num_branches ; i++) {
        mdpair a(initial+0,initial+length_of_branch*i+1);
        bindpairs.push_back(a);
    }

    for (int i = 0; i < num_branches; i++)
    {
        mdtriplet a(initial+0, initial+length_of_branch * i + 1, initial+length_of_branch*1+2);
        bendpairs.push_back(a);
    }

    for(int nb = 0 ; nb < num_branches ; nb++)
        for(int i = 1 ; i < length_of_branch ; i++ ) {
            mdpair a(initial + nb * length_of_branch + i, initial + nb * length_of_branch + i + 1);
            bindpairs.push_back(a);
        }

    for (int nb = 0; nb < num_branches; nb++)
        for (int i = 1; i < length_of_branch - 1; i++)
        {
            mdtriplet a(initial + nb * length_of_branch + i, initial + nb * length_of_branch + i + 1, initial + nb * length_of_branch + i + 2);
            bendpairs.push_back(a);
        }
    total_particles +=num_branches*length_of_branch+1;

}


matrix<int> convert_vector_list(const vector<mdpair> &x) {
    matrix<int> a(x.size(),2);
    #pragma omp parallel for
    for(int i = 0  ; i < x.size() ; i++) {
        a(i, 0) =  x[i].a;
        a(i, 1) = x[i].b;
    }
    return a;
}

void Nanostar::gets(matrix<int> &pairs, matrix<int> &specials, matrix<int> &not_specials){

    vector<mdpair> special_pairs;
    vector<mdpair> not_special_pairs;

    int totn  = pairs.getNsafe();

    #pragma omp parallel
    {
        vector<mdpair> index1_private;
        vector<mdpair> index2_private;

        index1_private.reserve(totn);
        index2_private.reserve(totn);
// vector<int> index3_private;
        #pragma omp for nowait schedule(static)
        for (int i = 0; i < pairs.getNsafe(); i++)
        {
            int p1 = pairs(i,0);
            int p2 = pairs(i,1);

            int which_nanostar = floor(p1 / total_parts_per_nanostar);
            int which_nanostar2 = floor(p2 / total_parts_per_nanostar);

            mdpair temp(p1, p2);

            if (which_nanostar != which_nanostar2 && count(stickerList.begin(), stickerList.end(), p1) > 0 && count(stickerList.begin(), stickerList.end(), p2) > 0)
            {
                //check if there are clashes
                // bool clash = false;
                // for(int j = 0  ; j < special_pairs.size() ; j++) {
                //     if(temp.either(special_pairs[j])) {
                //         mdpair temp2 = special_pairs[j];
                //         clash = true;
                //         int rr = rand() %2;
                //         if(rr == 0 ) {}
                //         else {special_pairs[j] = temp; }
                //     }
                // }
                // if(!clash)
                // {
                    index1_private.push_back(temp);
                // }
            }
            else{
               index2_private.push_back(temp);
            }

        }

        #pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); i++)
        {
        #pragma omp ordered
            special_pairs.insert(special_pairs.end(), index1_private.begin(), index1_private.end());
        }
        #pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); i++)
        {
        #pragma omp ordered
            not_special_pairs.insert(not_special_pairs.end(), index2_private.begin(), index2_private.end());
        }
    }



    specials = convert_vector_list(special_pairs);

    not_specials = convert_vector_list(not_special_pairs);

    //return pairs;
    // vector<mdpair> special_pairs;
    // vector<mdpair> not_special_pairs;
    // for(int i = 0 ; i < pairs.getNsafe() ; i++) {
    //     int p1 = pairs(i,0);
    //     int p2 = pairs(i,1);

    //     mdpair temp(p1,p2);

    //     int which_nanostar =  floor(p1/(num_branches*length_of_branch+1));
    //     int which_nanostar2 =  floor(p2/(num_branches*length_of_branch+1));

    //     bool endp1 = ((p1 - 1) % 10 == 0);
    //     bool endp2 = ((p2 - 1) % 10 == 0);

    //     if(which_nanostar != which_nanostar2  && endp1 && endp2) {
    //         special_pairs.push_back(temp);
    //     }
    //     else{
    //         not_special_pairs.push_back(temp);
    //     }


    // }
    }

void Nanostar::DoAnMC()
{

    int N = obj->getN();


    HSPotential wsa(10.0, 1.0);
    //HarmonicPotential spr(10., 1.1);
    FENEPotential spr(50, 1.4);

    matrix<int> bindones(N, 7);
    vector1<int> iterators(N);

    matrix<int> bp2(bindpairs.size(), 2);

    matrix<double> posi = obj->getdat();
    //matrix<double> centroids(num_nanostars,3); //find all the closest centroids


    
    // Collect all the interactions

    for (int i = 0; i < bindpairs.size(); i++)
    {
        mdpair temp = bindpairs[i];
        bp2(i, 0) = temp.a;
        bp2(i, 1) = temp.b;
    }

    for (int i = 0; i < bindpairs.size(); i++)
    {
        int p1 = bindpairs[i].a;
        int p2 = bindpairs[i].b;

        int iter1 = iterators[p1];
        int iter2 = iterators[p2];

        bindones(p1, iter1) = p2;
        bindones(p2, iter2) = p1;

        iterators[p1]++;
        iterators[p2]++;
    }
    cube newgeo = obj->getgeo();


    matrix<double> energy_per_bond(N, 7);
    for (int i = 0; i < N; i++)
    {
        int p1 = i;
        for (int j = 0; j < iterators[p1]; j++)
        {
            int p2 = bindones(i, j);

            double dis = newgeo.distance(posi, p1, p2);

            double en = wsa.energy(dis) + spr.energy(dis);

            energy_per_bond(i, j) = en;
        }
    }

    int ccc;
    int num = floor(l / 4.);

    matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num, ccc);

    int tb = boxes.getnrows();

    vector< vector<int> > each_box(tb);
   



    int cubes_per_length = (int)round(exp(log(tb) / 3));
    vector1<int> dim(dimension);
    for (int i = 0; i < dimension; i++)
    {
        int ij = 1;
        for (int j = 0; j < i; j++)
        {
            ij *= cubes_per_length;
        }
        dim[i] = ij;
    }

    for(int i = 0  ; i < N ; i++) {
        int c = newgeo.assign_box(posi,i,dim,cubes_per_length);

        each_box[c].push_back(i);
    }

    int ret = N/3;

    int No_Moves = 10000000;
    double step_size = 0.3;
    int downhill = 0;
    int uphill = 0;
    int reject = 0;

    vector1<double> energ(N);


        for (int i = 0; i < No_Moves; i++)
        {
            cout << i << endl;

            int choice = rand() % N; // choose a random particle

            double x = posi(choice, 0);
            double y = posi(choice, 1);
            double z = posi(choice, 2);

            //what box am I

            int wbami = newgeo.assign_box(posi,choice,dim,cubes_per_length);

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
            double aftE = 0;
            for (int j = 0; j < iterators[choice]; j++)
            {
                int p2 = bindones(choice, j);
                double dis1 = newgeo.distance(posi.getrowvector(choice), posi.getrowvector(p2));
                double dis2 = newgeo.distance(newpos, posi.getrowvector(p2));

                // double en = /* wsa.energy(dis) + */ spr.energy(dis);
                double enb = spr.energy(dis1);

                double ena = spr.energy(dis2);
                
                befE += enb;
                aftE += ena;
                // energy_per_bond(i, j) = en;
            }

            for(int box_no = 0 ; box_no < boxes.getncols() ; box_no++ ) {
                int p1 = choice;

                int boxunderc = boxes(wbami, box_no);
                
                for (int j = 0  ; j < each_box[boxunderc].size() ; j++) {
                    int p2 = each_box[boxes(wbami, box_no)][j];

                    if(p1 != p2) {
                    double dis1 = newgeo.distance(posi.getrowvector(p1), posi.getrowvector(p2));

                    double dis2 = newgeo.distance(newpos, posi.getrowvector(p2));

                    double enb = wsa.energy(dis1);

                    double ena = wsa.energy(dis2);


                    befE += enb;
                    aftE += ena;
                    }
                }
                      
            }

            


            double rf = (double)(rand()) / (double)(RAND_MAX);
            double dE = aftE - befE;


            double kT = 1.;
            if (dE < 0)
            {

                posi(choice, 0) = newx;
                posi(choice, 1) = newy;
                posi(choice, 2) = newz;
                // for (int j = 0; j < iterators[choice]; j++)
                // {
                //     energy_per_bond(choice, j) = es[j];
                // }
                downhill++;
                energ[choice] = aftE;

            }
            if (exp(-dE / kT) > rf)
            {
                // if (aftE > 1E10)
                // {
                //     cout << "uphill" << endl;
                //     cout << posi.getrowvector(choice) << endl;
                //     cout << newpos << endl;
                //     cout << choice << " " << befE << " " << aftE << endl;
                //     cout << "FENE" << endl;
                //     for (int j = 0; j < iterators[choice]; j++)
                //     {
                //         int p2 = bindones(choice, j);
                //         double dis1 = newgeo.distance(posi.getrowvector(choice), posi.getrowvector(p2));
                //         double dis2 = newgeo.distance(newpos, posi.getrowvector(p2));

                //         // double en = /* wsa.energy(dis) + */ spr.energy(dis);
                //         double enb = spr.energy(dis1);

                //         double ena = spr.energy(dis2);

                //         cout << p2 << " " << enb << " " << ena << endl;
                //         // energy_per_bond(i, j) = en;
                //     }
                //     cout << "LJ" << endl;
                //     for (int box_no = 0; box_no < boxes.getncols(); box_no++)
                //     {
                //         int p1 = choice;

                //         int boxunderc = boxes(wbami, box_no);

                //         for (int j = 0; j < each_box[boxunderc].size(); j++)
                //         {
                //             int p2 = each_box[boxes(wbami, box_no)][j];

                //             if (p1 != p2)
                //             {
                //                 double dis1 = newgeo.distance(posi.getrowvector(p1), posi.getrowvector(p2));

                //                 double dis2 = newgeo.distance(newpos, posi.getrowvector(p2));

                //                 double enb = wsa.energy(dis1);

                //                 double ena = wsa.energy(dis2);
                //                 cout << p2 << " " << enb << " " << ena << endl;
                //             }
                //         }
                //     }
                //     pausel();
                // }


                posi(choice, 0) = newx;
                posi(choice, 1) = newy;
                posi(choice, 2) = newz;

                // for (int j = 0; j < iterators[choice]; j++)
                // {
                //     energy_per_bond(choice, j) = es[j];
                // }
                uphill++;
                energ[choice] = aftE;
    
            }
            else
            {
                energ[choice] = befE;
                reject++;
            }

            if(i % ret == 0 && i > 0) {
                for(int j = 0 ; j < tb ; j++) {
                    each_box[j].clear(); //clear every box
                }
                each_box.clear();
                for(int j = 0 ; j < tb ; j++) {
                    each_box.push_back(vector<int>());
                }



                for (int j = 0; j < N; j++)
                {
                    int c = newgeo.assign_box(posi, j, dim, cubes_per_length); //repopulate every box

                    each_box[c].push_back(j);
                }
                
            }

            // double gettheta = atan2(sqrt(x ^ 2 + y ^ 2), z);
            // double getphi = atan2(y,x);

            // double newgettheta = gettheta + r1;
            // double newgetphi = gettheta + r2;
        }




        //try to relax the remaining high energy configurations
        step_size = 1.0;
        vector<int> high_energies;
        for(int i = 0  ; i < N ; i++) {
            if(energ[i] > 1E5) {
                high_energies.push_back(i);
            }
        }



        WCAPotential wsa2(10.0, 1.0,0.0);

        for (int i = 0; i < No_Moves; i++)
        {
            cout << i << endl;

            int choice2 = rand() % high_energies.size(); // choose a random particle
            int choice = high_energies[choice2];
            double x = posi(choice, 0);
            double y = posi(choice, 1);
            double z = posi(choice, 2);


            // what box am I

            int wbami = newgeo.assign_box(posi, choice, dim, cubes_per_length);

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
            double aftE = 0;
            for (int j = 0; j < iterators[choice]; j++)
            {
                int p2 = bindones(choice, j);
                double dis1 = newgeo.distance(posi.getrowvector(choice), posi.getrowvector(p2));
                double dis2 = newgeo.distance(newpos, posi.getrowvector(p2));

                // double en = /* wsa.energy(dis) + */ spr.energy(dis);
                double enb = spr.energy(dis1);

                double ena = spr.energy(dis2);

                befE += enb;
                aftE += ena;
                // energy_per_bond(i, j) = en;
            }



            for (int box_no = 0; box_no < boxes.getncols(); box_no++)
            {
                int p1 = choice;

                int boxunderc = boxes(wbami, box_no);

                for (int j = 0; j < each_box[boxunderc].size(); j++)
                {
                    int p2 = each_box[boxes(wbami, box_no)][j];

                    if (p1 != p2)
                    {
                        double dis1 = newgeo.distance(posi.getrowvector(p1), posi.getrowvector(p2));

                        double dis2 = newgeo.distance(newpos, posi.getrowvector(p2));

                        double enb = wsa2.energy(dis1);

                        double ena = wsa2.energy(dis2);

                        befE += enb;
                        aftE += ena;
                    }
                }
            }


            double rf = (double)(rand()) / (double)(RAND_MAX);
            double dE = aftE - befE;

            double kT = 1.;
            if (dE < 0)
            {
                posi(choice, 0) = newx;
                posi(choice, 1) = newy;
                posi(choice, 2) = newz;

                downhill++;
                energ[choice] = aftE;
            }
            if (exp(-dE / kT) > rf)
            {
                posi(choice, 0) = newx;
                posi(choice, 1) = newy;
                posi(choice, 2) = newz;

                uphill++;
                energ[choice] = aftE;
            }
            else
            {
                energ[choice] = befE;
                reject++;
            }

            if (i % ret == 0 && i > 0)
            {
                for (int j = 0; j < tb; j++)
                {
                    each_box[j].clear(); // clear every box
                }
                each_box.clear();
                for (int j = 0; j < tb; j++)
                {
                    each_box.push_back(vector<int>());
                }

                for (int j = 0; j < N; j++)
                {
                    int c = newgeo.assign_box(posi, j, dim, cubes_per_length); // repopulate every box

                    each_box[c].push_back(j);
                }
            }

            // double gettheta = atan2(sqrt(x ^ 2 + y ^ 2), z);
            // double getphi = atan2(y,x);

            // double newgettheta = gettheta + r1;
            // double newgetphi = gettheta + r2;
        }


        for (int choice = 0; choice < N; choice++)
        {
            int p1 = choice;
            int wbami = newgeo.assign_box(posi, choice, dim, cubes_per_length);
            for (int box_no = 0; box_no < boxes.getncols(); box_no++)
            {

                int boxunderc = boxes(wbami, box_no);

                for (int j = 0; j < each_box[boxunderc].size(); j++)
                {
                    int p2 = each_box[boxes(wbami, box_no)][j];

                    if (p1 != p2)
                    {
                        double dis1 = newgeo.distance(posi.getrowvector(p1), posi.getrowvector(p2));
                        if (dis1 < 0.8) {
                            cout << "(" << p1 << ": " << posi.getrowvector(p1) << ",  " << p2 <<": " << posi.getrowvector(p2) << "," << dis1 << ")"
                                 << ",";
                        }
                    }
                }
            }
        }
        //calculate the total energy



        obj->setdat(posi);

}

void Nanostar::run(int runtime, int every, int st = 0, string strbase = "")
{
    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);


    matrix<int> bp2(bindpairs.size(),2);

 
    matrix<int> bep2(bendpairs.size(), 3);

    //Collect all the interactions

    for(int i = 0 ; i < bindpairs.size() ; i++ ) {
        mdpair temp =  bindpairs[i];
        bp2(i, 0 ) = temp.a;
        bp2(i, 1) = temp.b;
    }

    for (int i = 0; i < bendpairs.size(); i++)
    {
        mdtriplet temp = bendpairs[i];
        bep2(i, 0) = temp.a;
        bep2(i, 1) = temp.b;
        bep2(i, 2) = temp.c;
    }


    int totalN = obj->getN();

    //int sibdiv = floor(ll/4.0);
    int ccc;
    int num = floor(l/4.);

    matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *froyo1 = obj->calculatepairs_parallel(boxes, 3.5);

    // *possible_stickers = new matrix<int>;
    matrix<int> possible_stickers;
    matrix<int> hs_pairs;




    this->gets(*froyo1,possible_stickers,hs_pairs);
    


     //matrix<int> *hs_pairs = new matrix<int>;
    //matrix<int> hs_pairs = this->get_non_s(*froyo1);

    unsigned int i;
     for (i = 0; i < runtime; i++)
     {
         cout << "r: " << i << endl;
         if (i > 0  && i % 25 == 0)
         {
             delete froyo1;

             // cout << "updated after: " << i << endl;
             // state = obj->getdat();
             froyo1 = obj->calculatepairs_parallel(boxes, 3.5);
             this->gets(*froyo1, possible_stickers, hs_pairs);
            //cout << "calculated" << endl;
         }

         

         matrix<double> F1((*obj).calculateforces(bp2, *bindp));

         //matrix<double> F2((*obj).calculateforces_threebody(bep2, *bendp));


         matrix<double> F3((*obj).calculateforces(possible_stickers, *faa));


         matrix<double> F4((*obj).calculateforces(hs_pairs, *hs));



        matrix<double> F = F1 /* + F2 + */+ F3 + F4;

        matrix<double> R(totalN, dimension);

        for (int i1 = 0; i1 < totalN; i1++)
        {
            for (int j = 0; j < dimension; j++)
            {
                R(i1, j) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
            }
         }
         if (i % every == 0)
         {

             //cout << i << endl;

             stringstream ss;

             ss << setw(number_of_digits) << setfill('0') << st+(i / every);

             matrix<double> pos = obj->getdat();

             string poss = "pos";
             poss = poss + strbase;

             poss += "_i=";

             string extension = ".csv";

             poss += ss.str();

             poss += extension;

             ofstream myfile;
             myfile.open(poss.c_str());

             myfile <<= pos;

             myfile.close();
         }

         (*obj).advance_mom(F, R);

         (*obj).advance_pos();
    }
}

#endif /* NANOSTAR_CPP */
