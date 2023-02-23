#ifndef PARTICLEADDER_H
#define PARTICLEADDER_H

#include "../MDBase/MD.h"

template <typename T> // requires c++11 compiler
T remove_at(std::vector<T> &v, typename std::vector<T>::size_type n)
{
    T ans = std::move_if_noexcept(v[n]);
    v[n] = std::move_if_noexcept(v.back());
    v.pop_back();
    return ans;
}

struct volume_vol {

    double ll;
    virtual bool in_volume(vector1<double> &myvec)=0;
    double get_ll() { return ll; }
    virtual volume_vol *clone() const = 0;
};

struct sphere_vol : public volume_vol  {
    double r = 5.; //default to 5,;
    double c1 =10.;
    double c2 = 10.;
    double c3 = 10.;

    bool in_volume(vector1<double> &myvec) {
        double tot = 0.0;

        tot += SQR(myvec[0] - c1);
        tot += SQR(myvec[1] - c2);
        tot += SQR(myvec[2] - c3);


        if(tot>SQR(r)) return false;
        else return true;

    }
    sphere_vol *clone() const
    {
        return new sphere_vol(*this);
    }
};


struct particle_adder {

    //this is a struct which has a list of indices, and weights in which they are added

    //we want it to return index i with weight [i], and then also to give a position to add the particle. 

    vector<double> weights;
    vector<int> indices;
    double rate; //this is the base rate at which particles are added, a particle is added for sure at this rate

    volume_vol *vol;

    particle_adder();

    void set_volume(volume_vol&);
    void set_weights(vector<double>&);
    void set_indices(vector<int>&);
    void set_rate(double &r);

    ~particle_adder() {delete vol;}

    vector1<double> generate_point_in_volume();

    template<class vec>
    void add_p(MD &obj, vec ind, bool&, vector1<double>&, int&); 
    //returns a particle at int with position vector1<> avoiding ind of dat




};

#include "ParticleAdder.cpp"

#endif /* NANOTUBE_H */