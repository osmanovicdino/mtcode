

Condensate::Condensate(double ll, int N, int M) : potential_bundle(vector1<potentialtheta3D*>(M)) {
    num = floor(ll/4.);

    obj = new LangevinNVTR;

    vector1<bool> pb(3,true);
    cube geo(ll,pb,3);

    matrix<double> dat(N,3);

    int pp = floor(ll-1.);

    vector<double> possible_pos_x;
    vector<double> possible_pos_y;
    vector<double> possible_pos_z;

    for(int i = 0 ; i < pp ; i++) {
        for(int j = 0  ; j < pp ; j++) {
            for(int k = 0 ; k < pp ; k++) {
                double x = 0.5 + i;
                double y = 0.5 + j;
                double z = 0.5 + k;
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

    double viscosity = 10.0;

    double hdradius = 0.5;

    b.setgamma(6.*pi*viscosity*hdradius);
    b.setgammar(8.*viscosity*hdradius*hdradius*hdradius);

    b.setkT(1.0);

    *obj = b;

}

void Condensate::set_potential_bundle(vector1<potentialtheta3D*> &a) {

int q =  potential_bundle.getsize();
for(int i = 0 ; i < q ; i++) {
    delete potential_bundle[i];
}
// cout << "old vals deleted" << endl;

potential_bundle = vector1<potentialtheta3D*>(a.getsize());

for(int i = 0 ; i < a.getsize() ; i++) {
    potential_bundle[i] = (a[i])->clone();
}

}

void Condensate::run(int runtime, int every, string strbase = "") {
    int ccc;

    int tf = ceil( (double)runtime / (double)every);
    int number_of_digits = 0;
    do {++number_of_digits;
    tf /= 10;
    } while(tf);



    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *pairs = obj->calculatepairs(boxes, 3.5);

    WCAPotential wsa(1.0,1.0,0.0);

    int NN = obj->getN();

    matrix<double> F(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> zeromatrix(NN,3);

    obj->calculateforces(*pairs, wsa);
    obj->calculate_forces_and_torques3D(*pairs, potential_bundle, F, T);

    obj->create_forces_and_torques_sphere(F, T);


    for (int i = 0; i < runtime; i++)
    {
        cout << i << endl;
        obj->advancemom_halfstep(F, T);
        obj->advance_pos();
        obj->rotate();

        F = obj->calculateforces(*pairs, wsa);
        
        T.reset(0.0);

        obj->calculate_forces_and_torques3D(*pairs, potential_bundle, F, T);

        obj->create_forces_and_torques_sphere(F, T);

        obj->advancemom_halfstep(F, T);

        if (i % every == 0)
        {
            delete pairs;
            pairs = obj->calculatepairs(boxes, 3.5);

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