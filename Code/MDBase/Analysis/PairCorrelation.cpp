#ifndef PAIRCORRELATION_CPP
#define PAIRCORRELATION_CPP

#include "../geometry.h"
#include "../../DataStructures/vector1.h"
#include "../../DataStructures/matrix2.h"


vector1<int> calcg1(matrix<double> &data, cube &geo, double max, double dx)
{
    vector<double> index1;

    double ll = geo.l;
    int ccc;
    int num = floor(ll / max); //NB THIS ALGORITHM GOES WRONG IF THERE ARE LESS THAN 3 BOXES num=2

    matrix<int> boxlist = geo.generate_boxes_relationships(num, ccc);

    int total_cubes = boxlist.getNsafe();

    double dims = (double)(geo.dimension);

    int cubes_per_length = (int)round(exp(log(total_cubes) / dims));

    vector<vector<int>> b;
    b.reserve(total_cubes);
    for (int j = 0; j < total_cubes; j++)
    {
        vector<int> temp;
        b.push_back(temp);
    }
    int dimension = (geo.dimension);

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

    for (int i = 0; i < data.getNsafe(); i++)
    {

        int c = geo.assign_box(data, i , dim, cubes_per_length);

        b[c].push_back(i);
    }

    int ss = boxlist.getncols();

#pragma omp parallel
    {
        vector<double> index1_private;

#pragma omp for nowait schedule(static)
        for (int c1 = 0; c1 < boxlist.getNsafe(); c1++)
            for (int c2 = 0; c2 < ss; c2++)
            {

                int box1 = c1;
                int box2 = boxlist(c1, c2);
                if (box1 == box2)
                {
                    for (int i = 0; i < b[box1].size(); i++)
                    {
                        for (int j = i + 1; j < b[box2].size(); j++)
                        {
                            int iterator1 = (b[box1])[i];
                            int iterator2 = (b[box2])[j];
                            //for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
                            //int q = ints.get_potential_number(iterator1,iterator2,l);

                            //if(int_dl[q]) {
                            double dis = geo.distance(data, iterator1, iterator2);

                            //int q = ints.get_potential_number(iterator1,iterator2,l);
                            index1_private.push_back(dis);
                            //index2_private.push_back(iterator2);
                            //index3_private.push_back(l);
                        }
                    }
                }
                else if (box2 > box1)
                {
                    for (int i = 0; i < b[box1].size(); i++)
                    {
                        for (int j = 0; j < b[box2].size(); j++)
                        {
                            int iterator1 = (b[box1])[i];
                            int iterator2 = (b[box2])[j];

                            double dis = geo.distance(data, iterator1, iterator2);

                            //int q = ints.get_potential_number(iterator1,iterator2,l);
                            index1_private.push_back(dis);
                            // for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
                            // 	int q = ints.get_potential_number(iterator1,iterator2,l);

                            // 	if(int_dl[q]) {
                            // 		bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);

                            // 		if(cond) {
                            // 			index1_private.push_back(iterator1);
                            // 			index2_private.push_back(iterator2);
                            // 			index3_private.push_back(l);
                            // 		}

                            // 	}
                            // 	else {
                            // 		index1_private.push_back(iterator1);
                            // 		index2_private.push_back(iterator2);
                            // 		index3_private.push_back(l);
                            // 	}
                            // }
                        }
                    }
                }
                else
                {
                }
            }

#pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); i++)
        {
#pragma omp ordered
            index1.insert(index1.end(), index1_private.begin(), index1_private.end());
        }
    }

    vector1<double> a(index1.size());
    //s_matrix<int> pairs(index1.size(),3);
    for (int i = 0; i < a.getsize(); i++)
    {
        a[i] = index1[i];
    }
    //cout << "new" << endl;
    //create  a histogram of all the distances
    int gs = floor(max / dx);
    vector1<int> g(gs);
    //where index i is (i)*dx to (i+1)*dx

    for (int i = 0; i < a.getsize(); i++)
    {
        int in = floor(a[i] / dx);
        if (in < gs)
            g[in] += 1;
    }

    return g;
}

double angd(vector1<double> &angs, int iterator1, int iterator2, double lim) {
    double asd1 = fmod(angs[iterator1], lim);
    double asd2 = fmod(angs[iterator2], lim);

    if(asd1 < 0)
        asd1 = 2 * pi + asd1;
    if (asd2 < 0)
        asd2 = 2 * pi + asd2;

    double dx = asd1 - asd2;

    if (abs(dx) > lim * 0.5)
        dx = dx - SIGN(lim, dx);

    // if(abs(dx) > lim/2.) { 
    //     cout << angs[iterator1] << endl;
    //     cout << angs[iterator2] << endl;
    //     cout << lim << endl;
    //     cout << asd1 << " " << asd2 << " " << dx << " " << iterator1 << " " << iterator2 << endl;
    //     error("asd");}
    return abs(dx);
}

matrix<int> calcg1(matrix<double> &data, vector1<double> &orient, cube &geo, double max, double dx, double dthe)
{
    vector<double> index1;
    vector<double> index2;

    double ll = geo.l;
    int ccc;
    int num = floor(ll / max); //NB THIS ALGORITHM GOES WRONG IF THERE ARE LESS THAN 3 BOXES num=2

    matrix<int> boxlist = geo.generate_boxes_relationships(num, ccc);

    int total_cubes = boxlist.getNsafe();

    double dims = (double)(geo.dimension);

    int cubes_per_length = (int)round(exp(log(total_cubes) / dims));

    vector<vector<int>> b;
    b.reserve(total_cubes);
    for (int j = 0; j < total_cubes; j++)
    {
        vector<int> temp;
        b.push_back(temp);
    }
    int dimension = (geo.dimension);

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

    for (int i = 0; i < data.getNsafe(); i++)
    {

        int c = geo.assign_box(data, i, dim, cubes_per_length);

        b[c].push_back(i);
    }

    int ss = boxlist.getncols();

    #pragma omp parallel
    {
        vector<double> index1_private;
        vector<double> index2_private;

        #pragma omp for nowait schedule(static)
        for (int c1 = 0; c1 < boxlist.getNsafe(); c1++)
            for (int c2 = 0; c2 < ss; c2++)
            {

                int box1 = c1;
                int box2 = boxlist(c1, c2);
                if (box1 == box2)
                {
                    for (int i = 0; i < b[box1].size(); i++)
                    {
                        for (int j = i + 1; j < b[box2].size(); j++)
                        {
                            int iterator1 = (b[box1])[i];
                            int iterator2 = (b[box2])[j];
                            //for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
                            //int q = ints.get_potential_number(iterator1,iterator2,l);

                            //if(int_dl[q]) {
                            double dis;// = geo.distance(data, iterator1, iterator2);
                            double dtheta = angd(orient, iterator1, iterator2, 2*pi);
                            vector1<double> un(2);
                            geo.distance_vector(data,iterator1,iterator2,un,dis);


                            double tempvx = cos(orient[iterator1]);
                            double tempvy = sin(orient[iterator1]);

                            vector1<double> as(2);
                            as[0] = tempvx;
                            as[1] = tempvy;

                            double val = scalar(as,un);

                            double angf = val/abs(val);
                            //int q = ints.get_potential_number(iterator1,iterator2,l);
                            index1_private.push_back(dis);
                            index2_private.push_back(angf*dtheta);
                            //index2_private.push_back(iterator2);
                            //index3_private.push_back(l);
                        }
                    }
                }
                else if (box2 > box1)
                {
                    for (int i = 0; i < b[box1].size(); i++)
                    {
                        for (int j = 0; j < b[box2].size(); j++)
                        {
                            int iterator1 = (b[box1])[i];
                            int iterator2 = (b[box2])[j];

                            double dis = geo.distance(data, iterator1, iterator2);
                            double dtheta = angd(orient, iterator1, iterator2, 2 * pi);
                            vector1<double> un(2);
                            geo.distance_vector(data, iterator1, iterator2, un, dis);

                            double tempvx = cos(orient[iterator1]);
                            double tempvy = sin(orient[iterator1]);

                            vector1<double> as(2);
                            as[0] = tempvx;
                            as[1] = tempvy;

                            double val = scalar(as, un);

                            double angf = val / abs(val);
                            //int q = ints.get_potential_number(iterator1,iterator2,l);
                            index1_private.push_back(dis);
                            index2_private.push_back(angf*dtheta);
                            // for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
                            // 	int q = ints.get_potential_number(iterator1,iterator2,l);

                            // 	if(int_dl[q]) {
                            // 		bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);

                            // 		if(cond) {
                            // 			index1_private.push_back(iterator1);
                            // 			index2_private.push_back(iterator2);
                            // 			index3_private.push_back(l);
                            // 		}

                            // 	}
                            // 	else {
                            // 		index1_private.push_back(iterator1);
                            // 		index2_private.push_back(iterator2);
                            // 		index3_private.push_back(l);
                            // 	}
                            // }
                        }
                    }
                }
                else
                {
                }
            }

#pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); i++)
        {
#pragma omp ordered
            index1.insert(index1.end(), index1_private.begin(), index1_private.end());
            index2.insert(index2.end(), index2_private.begin(), index2_private.end());
        }
    }

    vector1<double> a(index1.size());
    vector1<double> a2(index2.size());
    //s_matrix<int> pairs(index1.size(),3);
    for (int i = 0; i < a.getsize(); i++)
    {
        a[i] = index1[i];
        a2[i] = index2[i];
    }
    //create  a histogram of all the distances
    int gs = floor(max / dx);
    int ts = floor(2*pi / dthe);
    matrix<int> g(gs,ts);
    //where index i is (i)*dx to (i+1)*dx

    for (int i = 0; i < a.getsize(); i++)
    {
        int in = floor(a[i] / dx);
        int tn = floor((pi+a2[i])/dthe);
        if (in < gs) {
            // cout << endl;
            // cout << i << endl;
            // cout << in << " " << tn << endl;
            // cout << a2[i] << endl;
            g(in,tn) += 1;
        }
    }

    return g;
}
#endif /* PAIRCORRELATION_CPP */
