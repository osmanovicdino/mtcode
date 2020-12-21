#ifndef NANOSTAR_CPP
#define NANOSTAR_CPP

Nanostar::Nanostar(int N, double ll) :  bindpairs(vector<mdpair>()), bendpairs(vector<mdtriplet>()) {
    obj = new LangevinNVT;
    dimension = 3;
    l = ll;

    int num_nanostars = N;
    total_particles = 0;

    num_branches = 4;

    length_of_branch = 10;

    vector1<bool> is_periodic(dimension,true);
    cube bc(ll,is_periodic,dimension);

    double sigma = 1.0;

    double att_epp = 3.0;

    WCAPotential StickerPotential(10.0,sigma,att_epp);

    WCAPotential HS(10.0, sigma, 0.0);

    FENEPotential nrr(50,1.4);

    double preferred_angle = 0.0;

    double bending_strength = 10.0;

    BendingPotential nfr2(bending_strength,preferred_angle);

    for(int i =  0 ; i < num_nanostars; i++) {
        this->create_nanostar();
    }

    potential *q1 = StickerPotential.clone();
    faa = q1;
    potential *q2 = nrr.clone();
    bindp = q2;
    potential3 *q3 = nfr2.clone();
    bendp = q3;
    potential *q4 = HS.clone();
    hs = q4;


    matrix<double> store = create_initial_state();

    LangevinNVT b(bc);

    double kT = 1.0;
    double dt = 0.005;
    double eta = 50.;
    //gamma = eta;

    b.setdat(store);
    //b.setinteractions(nswca);
    b.setkT(kT);
    b.setdt(dt);
    b.setgamma(eta);
    b.setm(1.0);
    matrix<double> moms(b.getN(), dimension);
    b.setmom(moms);

    *obj = b;
}

matrix<double> Nanostar::create_initial_state() {

matrix<double> store(total_particles,3);

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

matrix<int> Nanostar::gets(matrix<int> &pairs, matrix<int> &specials, matrix<int> &not_specials){
    //return pairs;
    vector<mdpair> special_pairs;
    vector<mdpair> not_special_pairs;
    for(int i = 0 ; i < pairs.getNsafe() ; i++) {
        int p1 = pairs(i,0);
        int p2 = pairs(i,1);

        mdpair temp(p1,p2);

        int which_nanostar =  floor(p1/(num_branches*length_of_branch+1));
        int which_nanostar2 =  floor(p2/(num_branches*length_of_branch+1));

        bool endp1 = ((p1 - 1) % 10 == 0);
        bool endp2 = ((p2 - 1) % 10 == 0);

        if(which_nanostar != which_nanostar2  && endp1 && endp2) {
            special_pairs.push_back(temp);
        }
        else{
            not_special_pairs.push_back(temp);
        }


    }

}



void Nanostar::run(int runtime, int every, string strbase = "")
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
   
    matrix<int> *froyo1 = obj->calculatepairs(boxes, 3.5);

    // *possible_stickers = new matrix<int>;
    matrix<int> possible_stickers;
    matrix<int> hs_pairs;
    this->gets(*froyo1,possible_stickers,hs_pairs);

     //matrix<int> *hs_pairs = new matrix<int>;
    //matrix<int> hs_pairs = this->get_non_s(*froyo1);

    unsigned int i;
     for (i = 0; i < runtime; i++)
     {
         cout << i << endl;
         if (i % 25 == 0)
         {
             delete froyo1;

             // cout << "updated after: " << i << endl;
             // state = obj->getdat();
             froyo1 = obj->calculatepairs(boxes, 3.5);
             this->gets(*froyo1, possible_stickers, hs_pairs);
         }

         matrix<double> F1((*obj).calculateforces(bp2, *bindp));

         matrix<double> F2((*obj).calculateforces_threebody(bep2, *bendp));

         matrix<double> F3((*obj).calculateforces(possible_stickers, *faa));

         matrix<double> F4((*obj).calculateforces(hs_pairs, *hs));

         matrix<double> F = F1 + F2 + F3 + F4;

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

             ss << setw(number_of_digits) << setfill('0') << (i / every);

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
