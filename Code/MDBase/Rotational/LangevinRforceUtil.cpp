#ifndef LANGEVINRFORCEUTIL_CPP
#define LANGEVINRFORCEUTIL_CPP

matrix<int> LangevinNVTR::CreateEdgeList(matrix<int> &adj, vector1<int> &len)
{
    vector<mdpair> temp;
    temp.reserve(adj.getNsafe() * adj.getncols());

#pragma omp parallel
    {
        vector<mdpair> tempprivate;
        tempprivate.reserve(adj.getNsafe() * adj.getncols());
#pragma omp for nowait schedule(static)
        for (int i = 0; i < adj.getNsafe(); i++)
        {
            for (int j = 0; j < len[i]; j++)
            {
                int p1 = i;
                int p2 = adj(i, j);
                if (p2 > p1)
                {
                    mdpair mypair(p1, p2);
                    tempprivate.push_back(mypair);
                }
            }
        }

#pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); i++)
        {
#pragma omp ordered
            temp.insert(temp.end(), tempprivate.begin(), tempprivate.end());
        }
    }
    matrix<int> a(temp.size() * 2, 2);
//s_matrix<int> pairs(index1.size(),3);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < temp.size(); i++)
    {
        (a)(2 * i, 0) = temp[i].b;
        (a)(2 * i, 1) = temp[i].a;

        (a)(2 * i + 1, 0) = temp[i].a;
        (a)(2 * i + 1, 1) = temp[i].b;
    }
    return a;
}

void get_energies(int i1, int i2, int i3, int nb1, int nb2, int nb3, const matrix<int> &boindices, const matrix<double> &boenergies, double &e12, double &e23, double &e13) {

    for(int i = 0 ; i < nb1 ; i++) {
    if(boindices.gpcons(i1,i) == i2) {
        e12 = boenergies.gpcons(i1,i);
    }
    else if(boindices.gpcons(i1,i) == i3) {
        e13 = boenergies.gpcons(i1,i);
    }
    else{

    }

    }

    for (int i = 0; i < nb2; i++)
    {
        if (boindices.gpcons(i2, i) == i1)
        {
            e12 = boenergies.gpcons(i2, i);
        }
        else if (boindices.gpcons(i2, i) == i3)
        {
            e23 = boenergies.gpcons(i2, i);
        }
        else
        {
        }
    }

    for (int i = 0; i < nb3; i++)
    {
        if (boindices.gpcons(i3, i) == i1)
        {
            e13 = boenergies.gpcons(i3, i);
        }
        else if (boindices.gpcons(i3, i) == i2)
        {
            e23 = boenergies.gpcons(i3, i);
        }
        else
        {
        }
    }

    
}



vector<mdpairwd> RecapitulateEdges(int size_of_cluster, const matrix<int> &boindices2, const matrix<mdpairwd> &boindices, const vector1<int> &tempbound, int i, vector<int> &unique_indexes)
{
    vector<mdpairwd> matched;
    matched.reserve(size_of_cluster * (size_of_cluster - 1) / 2);

    
    for (int j = 0; j < size_of_cluster; j++)
    {
        int myindex = boindices2.gpcons(i, j);
        //unique_indices[j] = myindex;
        //for (int k = 0; k < tempbound[indexes[j]]; k++)
        for (int k = 0; k < tempbound.gpcons(myindex); k++)
        {
            mdpairwd m1;
            //int f1 = indexes[j];
            //int f2 = boindices(indexes[j], k);

            int f1 = myindex;
            int f2 = boindices.gpcons(myindex, k).a;

            //cout << f1 << " " << f2 << endl;

            std::vector<int>::iterator itr1 = std::find(unique_indexes.begin(), unique_indexes.end(), f1);
            std::vector<int>::iterator itr2 = std::find(unique_indexes.begin(), unique_indexes.end(), f2);

            //cout << std::distance(unique_indexes.begin(), itr1) << " " << std::distance(unique_indexes.begin(), itr2) << endl;
            f1 = std::distance(unique_indexes.begin(), itr1);
            f2 = std::distance(unique_indexes.begin(), itr2);

            if (f1 < f2)
            {
                m1.a = f1;
                m1.b = f2;
            }
            else
            {
                m1.a = f2;
                m1.b = f1;
            }

            m1.scr = boindices.gpcons(myindex, k).scr;


            matched.push_back(m1);
        }
    }
    sort(matched.begin(), matched.end());
    matched.erase(unique(matched.begin(), matched.end()), matched.end());

    return matched;
}

void RecapitulateEdges(int size_of_cluster, const matrix<int> &boindices2, const matrix<int> &boindices, const matrix<double> &boscores, const vector1<int> &tempbound, int i, vector<int> &unique_indexes, vector<mdpairwd> &matched)
{
    // vector<mdpairwd> matched;
    // matched.reserve(size_of_cluster * (size_of_cluster - 1) / 2);

    for (int j = 0; j < size_of_cluster; j++)
    {
        int myindex = boindices2.gpcons(i, j);
        //unique_indices[j] = myindex;
        //for (int k = 0; k < tempbound[indexes[j]]; k++)
        for (int k = 0; k < tempbound.gpcons(myindex); k++)
        {
            
            mdpairwd m1;
            //int f1 = indexes[j];
            //int f2 = boindices(indexes[j], k);

            int f1 = myindex;
            int f2 = boindices.gpcons(myindex, k);


            //cout << f1 << " " << f2 << endl;

            std::vector<int>::iterator itr1 = std::find(unique_indexes.begin(), unique_indexes.end(), f1);
            std::vector<int>::iterator itr2 = std::find(unique_indexes.begin(), unique_indexes.end(), f2);

            //cout << std::distance(unique_indexes.begin(), itr1) << " " << std::distance(unique_indexes.begin(), itr2) << endl;
            f1 = std::distance(unique_indexes.begin(), itr1);
            f2 = std::distance(unique_indexes.begin(), itr2);

            if (f1 < f2)
            {
                m1.a = f1;
                m1.b = f2;
            }
            else
            {
                m1.a = f2;
                m1.b = f1;
            }

            m1.scr = boscores.gpcons(myindex, k);

            int totmatched = matched.size();
           
            bool already = false;
            for(int lj = 0 ; lj < totmatched ; lj++) {
                int g1,g2;
                matched[lj].get(g1,g2);
                
                if (g1 == m1.a && g2 == m1.b) {
                    //cout << "already" << endl;
                    already = true;
                    break;
                }
            }


            if(!already) {
            //cout << "pushed back" << endl;
            matched.push_back(m1);
            }
        }
    }

    // sort(matched.begin(), matched.end());
    // matched.erase(unique(matched.begin(), matched.end()), matched.end());

}

void Process_Doublet(int ti1, int ti2, BinaryBindStore &bo, AbstractBindingModel &bm, vector<mdpair> &mypairs_private)
{
    int i1;
    int i2;
    sort_doublet(ti1, ti2, i1, i2);

    bool alreadybound_to_eachother = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];

    bool aft;

    bm.doublet(alreadybound_to_eachother, i1, i2, aft);

    if (aft)
    {
        bo.boundto[i1] = i2;
        bo.boundto[i2] = i1;
        bo.isbound[i1] = true;
        bo.isbound[i2] = true;
        //tbt++;
        mypairs_private.push_back(mdpair(i1, i2));
    }
    else
    {

        bo.isbound[i1] = false;
        bo.isbound[i2] = false;
    }
}

void Process_Triplet(int ti1, int ti2, int ti3, BinaryBindStore &bo, AbstractBindingModel &bm, matrix<mdpair> &boindices, vector1<int> &tempbound, vector<mdpair> &mypairs_private ) {
    int i1;
    int i2;
    int i3;

    sort_triplet(ti1, ti2, ti3, i1, i2, i3);

    //DETERMINE WHETHER THEY ARE BOUND
    bool b12 = bo.boundto[i1] == i2 && bo.boundto[i2] == i1 && bo.isbound[i1] && bo.isbound[i2];
    bool b23 = bo.boundto[i2] == i3 && bo.boundto[i3] == i2 && bo.isbound[i2] && bo.isbound[i3];
    bool b13 = bo.boundto[i1] == i3 && bo.boundto[i3] == i1 && bo.isbound[i1] && bo.isbound[i3];

    //DETERMINE THE CONNECTIVENESS OF THE GRAPH
    //remember, that in order to count as a triplet
    bool c12 = false;
    bool c23 = false;
    bool c13 = false;

    int nb1 = tempbound[i1];
    int nb2 = tempbound[i2];
    int nb3 = tempbound[i3];

    // if(nb1 > 2 || nb2 > 2 || nb3 > 2) {

    //     cout << i1 << " " << i2 << " " << i3 << endl;

    //     outfunc(tempbound, "t1");
    //     outfunc(boindices, "t2");

    //     pausel();
    //     error("something weird in code");
    // }

    if (nb1 == 1)
    {
        int tempi = boindices(i1, 0).a;
        if (tempi == i2)
        {
            c12 = true;
            c13 = false;
            c23 = true; //in order to be a triplet
        }
        else if (tempi == i3)
        {
            c13 = true;
            c12 = false;
            c23 = true;
        }
        else
            error("something weird");

        //check the other
    }
    else if (nb1 == 2)
    {
        c12 = true;
        c13 = true;

        int nb2 = tempbound[i2];
        if (nb2 == 1)
        {
            c23 = false;
        }
        else
        {
            c23 = true;
        }
    }
    else
    {
        //cout << ccs[i] << " " << endl;
        cout << 3 << endl;
        cout << b12 << " " << b23 << " " << b13 << endl;
        cout << i1 << " " << i2 << " " << i3 << endl;
        cout << tempbound[i1] << " " << tempbound[i2] << " " << tempbound[i3] << endl;

        //cout << indexes2[i1] << " " << indexes2[i2] << " " << indexes2[i3] << endl;

        // for (int k = 0; k < total_number_of_patches; k++)
        // {
        //     if (indexes2[k] == indexes2[i1])
        //     {
        //         cout << k << " ";
        //     }
        // }
        cout << endl;
        cout << "indexes done" << endl;
        cout << "parallel vs serial" << endl;
        cout << endl;

        error("error in triplet clustering algorithm");
    }

    bool a12;
    bool a23;
    bool a13;
    //cout << "triplet called: " << i1 << " " << i2 << " " << i3 << endl;
    // //pausel();
    // bool cond1 = i1 < 12000 && i2 > 12000 + 4 * 15000 && (i2 - 12000 - 15000 * 4) % 3 == 2;
    // bool cond2 = i1 < 12000 && i3 > 12000 + 4 * 15000 && (i3 - 12000 - 15000 * 4) % 3 == 2;
    // bool cond3 = i2 < 12000 && i3 > 12000 + 4 * 15000 && (i3 - 12000 - 15000 * 4) % 3 == 2;

    bm.triplet(b12, b23, b13, c12, c23, c13, i1, i2, i3, a12, a23, a13);

    if (a12)
    {
        bo.boundto[i1] = i2;
        bo.boundto[i2] = i1;
        bo.isbound[i1] = true;
        bo.isbound[i2] = true;
        bo.isbound[i3] = false;
        mypairs_private.push_back(mdpair(i2, i1));
    }
    else if (a23)
    {

        bo.boundto[i2] = i3;
        bo.boundto[i3] = i2;
        bo.isbound[i1] = false;
        bo.isbound[i2] = true;
        bo.isbound[i3] = true;
        mypairs_private.push_back(mdpair(i3, i2));
    }
    else if (a13)
    {
        bo.boundto[i1] = i3;
        bo.boundto[i3] = i1;
        bo.isbound[i1] = true;
        bo.isbound[i2] = false;
        bo.isbound[i3] = true;
        mypairs_private.push_back(mdpair(i3, i1));
    }
    else
    {
        bo.isbound[i1] = false;
        bo.isbound[i2] = false;
        bo.isbound[i3] = false;
    }
}

#endif /* LANGEVINRFORCEUTIL_CPP */
