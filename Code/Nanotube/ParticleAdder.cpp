particle_adder::particle_adder() : weights(vector<double>()),indices(vector<int>()) {
    rate = 0.0;
    vol = new sphere_vol;
}

void particle_adder::set_volume(volume_vol &v)
{
    delete vol;

    vol = v.clone();
}

void particle_adder::set_rate(double &r) {
    rate = r;
}

void particle_adder::set_weights(vector<double> &we)
{
    weights = we;
}

void particle_adder::set_indices(vector<int> &in)
{
    indices = in;
}

vector1<double> particle_adder::generate_point_in_volume() {
    double ll =  vol->get_ll();
    
    double r1 = ll * ((double)rand() / (double)(RAND_MAX));
    double r2 = ll * ((double)rand() / (double)(RAND_MAX));
    double r3 = ll * ((double)rand() / (double)(RAND_MAX));

    vector1<double> myvec1(3);
    myvec1[0] = r1;
    myvec1[1] = r2;
    myvec1[2] = r3;
    int iter = 0;
    while (!(vol->in_volume(myvec1)))
    {
        double r1 = ll * ((double)rand() / (double)(RAND_MAX));
        double r2 = ll * ((double)rand() / (double)(RAND_MAX));
        double r3 = ll * ((double)rand() / (double)(RAND_MAX));
        myvec1[0] = r1;
        myvec1[1] = r2;
        myvec1[2] = r3;
        iter++;
        //if(iter > 100) error("can't add new particle, the chosen volume is too small for this to be an efficient algorithm");
        
    }

    return myvec1;
}

template <class vec> 
void particle_adder::add_p(MD &obj, vec ind, bool &does_return, vector1<double> &ve, int &fi) {

does_return = false;
if(weights.size() <1 ) {

}
else {
double r1 = ((double)rand() / (double)(RAND_MAX));
    if(r1<rate) {
        
        does_return = true;
        vector1<double> myvec1(generate_point_in_volume());



        vector1<double> mydistances(ind.size());

        for (int j = 0; j < ind.size(); j++)
        {
            int i  =  ind[j];

            double x2 = obj.getcoordinate(i, 0);
            double y2 = obj.getcoordinate(i, 1);
            double z2 = obj.getcoordinate(i, 2);
            vector1<double> myvec2(3);
            myvec2[0] = x2;
            myvec2[1] = y2;
            myvec2[2] = z2;
            mydistances[j] = obj.getgeo().distance(myvec1, myvec2);
        }

        double mymin;
        if( ind.size() > 0) {
            mymin = minval(mydistances);
        }
        else {
            mymin=2.;
        }
        int iter = 0;
        
        while (mymin < 1.0)
        {
            iter++;
            myvec1 = generate_point_in_volume();
            
            for (int j = 0; j < ind.size(); j++)
            {
                int i = ind[j];

                double x2 = obj.getcoordinate(i, 0);
                double y2 = obj.getcoordinate(i, 1);
                double z2 = obj.getcoordinate(i, 2);
                vector1<double> myvec2(3);
                myvec2[0] = x2;
                myvec2[1] = y2;
                myvec2[2] = z2;
                mydistances[j] = obj.getgeo().distance(myvec1, myvec2);
            }
            mymin = minval(mydistances);

            if (iter > 1000) {
                cout << "couldn't place" << endl;
                pausel();
            }
        }

        // cout << mymin << endl;

        double result = std::reduce(weights.begin(), weights.end());

        double randv = result*((double)rand() / (double)(RAND_MAX));
        double w = 0.0;
        int index_which =0;
        for(int i = 0 ; i < indices.size() ; i++) {
        w += weights[i];
        if(w>randv) {
            index_which = i;
            break;
        }    
        }

  

        ve = myvec1;
        fi = indices[index_which];



        remove_at(weights,index_which);
        remove_at(indices, index_which);



    }
    else{
        //do nothing
    }
}

} 