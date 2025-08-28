#ifndef LANGEVINRFORCEUTIL_CPP
#define LANGEVINRFORCEUTIL_CPP

matrix<int> LangevinNVTR::CreateEdgeList(matrix<int> &adj, vector1<int> &len)
{
    vector<mdpair> temp;
    temp.reserve(adj.getNsafe() * adj.getncols());
    int tota = adj.getNsafe();
    #pragma omp parallel
    {
        vector<mdpair> tempprivate;
        tempprivate.reserve(adj.getNsafe() * adj.getncols());
        #pragma omp for nowait schedule(static)
        for (int i = 0; i < tota; i++)
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
    int tempn = temp.size();
    matrix<int> a(tempn * 2, 2);
//s_matrix<int> pairs(index1.size(),3);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < tempn; i++)
    {
        (a)(2 * i, 0) = temp[i].b;
        (a)(2 * i, 1) = temp[i].a;

        (a)(2 * i + 1, 0) = temp[i].a;
        (a)(2 * i + 1, 1) = temp[i].b;
    }
    return a;
}

template <class Q>
matrix<int> Bond_Count(BinaryBindStore &bbs, Q &sorter, int no_types) {
    int N = bbs.boundto.getsize();
    
    // vector1<int> counts((no_types)*(no_types-1)/2);

    matrix<int> counts(no_types,no_types);

    for(int i = 0 ; i < N ; i++) {
        if(bbs.isbound[i]) {
            int p1 =  sorter(i);
            int p2 =  sorter(bbs.boundto[i]);

            counts(p1-1,p2-1)++;
        }
    }

    return counts;
}

template <class Q>
matrix<int> Bond_Change(BinaryBindStore &bbs1, BinaryBindStore &bbs2, Q &sorter, int no_types) {
    int N = bbs1.boundto.getsize();
    //what bonds are changing, where bbs2 is the later time point
    matrix<int> counts( (2*no_types) , (no_types));
    //int neg = ((no_types) * (no_types + 1))/2;
    int neg = no_types;
    for (int i = 0; i < N; i++)
    {
        int i1;
        if( bbs1.isbound[i] == false && bbs2.isbound[i] == true ) {
            int p1 = sorter(i)-1;
            int p2 = sorter(bbs2.boundto[i])-1;

            counts(p1 , p2)++; //p1,p2 bond formed

        }
        else if (bbs1.isbound[i] == true && bbs2.isbound[i] == false )
        {
            int p1 = sorter(i)-1;
            int p2 = sorter(bbs1.boundto[i])-1;

            counts(neg + p1 , p2)++; //p1,p2 bond lost
        }
        else if (bbs1.isbound[i] == true && bbs2.isbound[i] == true && (bbs1.boundto[i] != bbs2.boundto[i]) )
        {
            int p1 = sorter(i)-1;
            int p2 = sorter(bbs1.boundto[i])-1;

            int p3 = sorter(i)-1;
            int p4 = sorter(bbs2.boundto[i])-1;

            counts(neg + p1 , p2)++; //p1,p2 bond lost
            counts(p3 , p4)++; //p3,p4 bond formed
            
        }
        else{
            //do nothing
        }


        // 0 - > 1
        // 1 -> 0
        // 1 -> 1 but i->j

        //what are the case for difference
        

    }

    return counts;
}

template <class Q>
void  Bond_Change_Type(BinaryBindStore &bbs1, BinaryBindStore &bbs2, matrix<int> &boindices2, vector1<int> &ccs, Q &my_sorter, int no_types) {

//go through the whole thing

matrix<int> doublets_formed(no_types,no_types);

matrix<int> doublets_lost(no_types, no_types);

matrix<int> triplets_formed(no_types, no_types);

matrix<int> triplets_lost(no_types, no_types);

matrix<int> nlets_formed(no_types, no_types);

matrix<int> nlets_lost(no_types, no_types);

for(int i = 0  ; i < ccs.getsize() ; i++) {

    if(ccs[i] == 2) {
        //were they bound before?

        int wp1 = boindices2(i,0);
        int wp2 =  boindices2(i,1);

        int t1 = my_sorter(wp1) - 1;
        int t2 = my_sorter(wp2) - 1;

        bool cond1 = bbs1.boundto[wp1] == wp2 && bbs1.boundto[wp2] == wp1;
        bool b1, b2;

        b1 = bbs1.isbound[wp1];
        b2 = bbs1.isbound[wp2];

        bool alreadybound =  b1 && b2 && cond1;

        bool cond2 = bbs2.boundto[wp1] == wp2 && bbs2.boundto[wp2] == wp1;
        bool b3, b4;

        b3 = bbs2.isbound[wp1];
        b4 = bbs2.isbound[wp2];

        bool afterbound =  b3 && b4 && cond2;

        if(alreadybound && afterbound) {

        }
        else if(alreadybound && !afterbound) {
            doublets_lost(t1,t2)++;
        }
        else if(!alreadybound && afterbound){
            doublets_formed(t1,t2)++;
        }
        else{

        }


    }

    else if(ccs[i] == 3 ) {

        int size_of_cluster = ccs[i];

        vector<mdpair> matched;
        matched.reserve(size_of_cluster * (size_of_cluster - 1) / 2);

        for (int j = 0; j < size_of_cluster; j++)
        {
            int myindex = boindices2(i, j);
            //for (int k = 0; k < tempbound[indexes[j]]; k++)
            for (int k = j + 1; k < size_of_cluster; k++)
            {
                mdpair m1;
                //int f1 = indexes[j];
                //int f2 = boindices(indexes[j], k);

                int f1 = myindex;
                int f2 = boindices2(i, k);

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
                matched.push_back(m1);
            }
        }


        vector1<bool> bounded_before(matched.size());

        for (int j = 0; j < matched.size(); j++)
        {
            int i1 = matched[j].a;
            int i2 = matched[j].b;
            bool b12 = bbs1.boundto[i1] == i2 && bbs1.boundto[i2] == i1 && bbs1.isbound[i1] && bbs1.isbound[i2];
            bounded_before[j] = b12;
        }

        vector1<bool> bounded_after(matched.size());

        for (int j = 0; j < matched.size(); j++)
        {
            int i1 = matched[j].a;
            int i2 = matched[j].b;
            bool b12 = bbs2.boundto[i1] == i2 && bbs2.boundto[i2] == i1 && bbs2.isbound[i1] && bbs2.isbound[i2];
            bounded_after[j] = b12;
        }

        for (int j = 0; j < matched.size(); j++)
        {
            int typea = my_sorter(matched[j].a) - 1;
            int typeb = my_sorter(matched[j].b) - 1;

            if (bounded_after[j] && !bounded_before[j])
            {
                triplets_formed(typea, typeb)++;
            }
            else if (!bounded_after[j] && bounded_before[j])
            {
                triplets_lost(typea, typeb)++;
            }
            else
            {
            }
        }

    }
    else if(ccs[i] > 3) {

        int size_of_cluster = ccs[i];

        vector<mdpair> matched;
        matched.reserve(size_of_cluster * (size_of_cluster - 1) / 2);

        for (int j = 0; j < size_of_cluster; j++)
        {
            int myindex = boindices2(i, j);
            //for (int k = 0; k < tempbound[indexes[j]]; k++)
            for (int k = j+1; k < size_of_cluster; k++)
            {
                mdpair m1;
                //int f1 = indexes[j];
                //int f2 = boindices(indexes[j], k);

                int f1 = myindex;
                int f2 = boindices2(i, k);

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
                matched.push_back(m1);
            }
        }

        vector1<bool> bounded_before(matched.size());

        for (int j = 0; j < matched.size(); j++)
        {
            int i1 = matched[j].a;
            int i2 = matched[j].b;
            bool b12 = bbs1.boundto[i1] == i2 && bbs1.boundto[i2] == i1 && bbs1.isbound[i1] && bbs1.isbound[i2];
            bounded_before[j] = b12;
        }

        vector1<bool> bounded_after(matched.size());

        for (int j = 0; j < matched.size(); j++)
        {
            int i1 = matched[j].a;
            int i2 = matched[j].b;
            bool b12 = bbs2.boundto[i1] == i2 && bbs2.boundto[i2] == i1 && bbs2.isbound[i1] && bbs2.isbound[i2];
            bounded_after[j] = b12;
        }

        for(int j = 0  ; j < matched.size() ; j++) {
            int typea = my_sorter(matched[j].a)-1;
            int typeb = my_sorter(matched[j].b)-1;

            if(bounded_after[j] && !bounded_before[j]) {
                nlets_formed(typea,typeb)++;
            }
            else if (!bounded_after[j] && bounded_before[j])
            {
                nlets_lost(typea , typeb )++;
            }
            else {

            }

        }


    }
    else{

    }

}

cout <<= doublets_formed;
cout <<= triplets_formed;
cout <<= nlets_formed;

cout <<= doublets_lost;
cout <<= triplets_lost;
cout <<= nlets_lost;


}

template <class Q>
bool is_there_a_13( int i1, int i2, int i3, Q &sorter)
{

    if ( (((sorter(i1) == 1 && sorter(i2) == 3) || (sorter(i1) == 3 && sorter(i2) == 1))))
    {
        return true;
    }
    else if ( ((sorter(i3) == 1 && sorter(i2) == 3) || (sorter(i2) == 3 && sorter(i3) == 1)))
    {
        return true;
    }
    else if ( ((sorter(i3) == 1 && sorter(i1) == 3) || (sorter(i1) == 3 && sorter(i3) == 1)))
    {
        return true;
    }
    else
    {
        return false;
    }
}


template <class Q>
bool is_there_a_13_bond(bool b12, bool b23, bool b13, int i1, int i2 , int i3, Q &sorter) {
     
    if(b12 && (( (sorter(i1)==1 && sorter(i2)==3) || (sorter(i1)==3 && sorter(i2)==1) )) ) {
        return true;
    }
    else if(b23 && ( (sorter(i3)==1 && sorter(i2)==3) || (sorter(i2)==3 && sorter(i3)==1) )) {
        return true;
    }
    else if(b13 && ( (sorter(i3)==1 && sorter(i1)==3) || (sorter(i1)==3 && sorter(i3)==1) ) ) {
        return true;
    }
    else {
        return false;
    }
}

template <class Q>
bool is_123(int i1, int i2, int i3, Q &sorter)
{

    int ti1 = sorter(i1);
    int ti2 = sorter(i2);
    int ti3 = sorter(i3);

    int k1,k2,k3;
    sort_triplet(ti1,ti2,ti3,k1,k2,k3);

    if(k1 == 1 && k2 == 2 && k3 ==3 ) {
        return true;
    }
    else {
        return false;
    }
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
    double r = (double)rand()/(double)RAND_MAX; //NEED TO FIX BEFORE USE

    bm.doublet(alreadybound_to_eachother, i1, i2, aft,r);

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
    double r = (double)rand()/(double)(RAND_MAX); //NEED TO FIX THIS IF I EVER USE THIS FUNCTION

    bm.triplet(b12, b23, b13, c12, c23, c13, i1, i2, i3, a12, a23, a13, r);

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

void LangevinNVTR::setup_random_binding(vector<patchint> &pairs, vector<int> &divs, ComboPatch &iny, BinaryBindStore &bo, AbstractBindingModel &bm)
{
    int total_number_of_patches = bo.boundto.getsize(); //iny.get_total_patches(this->getN());

    vector1<int> tempbound;
    tempbound.resize_parallel(total_number_of_patches); //no binding to begin wtih

    int depth_of_matrix = 10; //Choose this value to be deep enough such that all values can be stored

    matrix<int> boindices;   //(total_number_of_patches, depth_of_matrix);
    matrix<double> boscores; //(total_number_of_patches, depth_of_matrix);

    boindices.resize_parallel(total_number_of_patches, depth_of_matrix);
    boscores.resize_parallel(total_number_of_patches, depth_of_matrix);

    vector<mdpairwd> edgelist;
    edgelist.reserve(total_number_of_patches);

    //std::mutex mtx;

    //int total_checks = 0;
    int tn_pairs = pairs.size();
    int t_u_pairs = divs.size();
    double mk = iny.max_check;

    if (tn_pairs > 0)
    {
#pragma omp parallel
        {
            vector<mdpairwd> edgelist_private;
            edgelist_private.reserve(total_number_of_patches);

#pragma omp for nowait schedule(dynamic)
            for (int ik = 0; ik < t_u_pairs + 1; ++ik)
            {
                int i, fi; //the start and end indices
                if (ik == 0)
                {
                    if (tn_pairs == 0)
                    {
                        i = 0;
                        fi = 0;
                        error("should get here");
                    }
                    else if (t_u_pairs == 0)
                    {
                        i = 0;
                        fi = tn_pairs;
                    }
                    else
                    {
                        i = 0;
                        fi = divs[ik];
                    }
                }
                else if (ik == t_u_pairs)
                {
                    i = divs[ik - 1];
                    fi = tn_pairs;
                }
                else
                {
                    i = divs[ik - 1];
                    fi = divs[ik];
                }

                int p1, p2;
                pairs[i].get_particle(p1, p2);
                // patchint tempo = pairs[i];
                // int p1 = tempo.particle_index1;
                // int p2 = tempo.particle_index2;

                //int i1 = pairs(i,2);
                double dis;
                //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
                vector1<double> un(dimension);
                geo.distance_vector(*dat, p1, p2, un, dis);

                //un = i-j

                if (dis < SQR(mk))
                {
                    dis = sqrt(dis);
                    un /= dis;
                    double dx = un.gpcons(0);
                    double dy = un.gpcons(1);
                    double dz = un.gpcons(2);

                    double qtemp0 = orient->gpcons(p1, 0);
                    double qtemp1 = orient->gpcons(p1, 1);
                    double qtemp2 = orient->gpcons(p1, 2);
                    double qtemp3 = orient->gpcons(p1, 3);
                    double qtemp4 = orient->gpcons(p1, 4);
                    double qtemp5 = orient->gpcons(p1, 5);
                    double qtemp6 = orient->gpcons(p1, 6);
                    double qtemp7 = orient->gpcons(p1, 7);
                    double qtemp8 = orient->gpcons(p1, 8);

                    double gtemp0 = orient->gpcons(p2, 0);
                    double gtemp1 = orient->gpcons(p2, 1);
                    double gtemp2 = orient->gpcons(p2, 2);
                    double gtemp3 = orient->gpcons(p2, 3);
                    double gtemp4 = orient->gpcons(p2, 4);
                    double gtemp5 = orient->gpcons(p2, 5);
                    double gtemp6 = orient->gpcons(p2, 6);
                    double gtemp7 = orient->gpcons(p2, 7);
                    double gtemp8 = orient->gpcons(p2, 8);

                    for (int j = i; j < fi; j++)
                    {
                        //pairschecked++;
                        int potn, wp1, wp2;
                        pairs[j].get_patch(potn, wp1, wp2);
                        mypot *temppot = iny.potential_bundle[potn];

                        double nxb1 = temppot->nxb1;
                        double nxb2 = temppot->nxb2;
                        double nyb1 = temppot->nyb1;
                        double nyb2 = temppot->nyb2;
                        double nzb1 = temppot->nzb1;
                        double nzb2 = temppot->nzb2;
                        double disp = temppot->interaction_distance;
                        double thetam = temppot->thetam;
                        double ctm = cos(thetam);
                        double nx1 = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
                        double ny1 = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
                        double nz1 = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;

                        double nx2 = nxb2 * gtemp0 + nyb2 * gtemp3 + nzb2 * gtemp6;
                        double ny2 = nxb2 * gtemp1 + nyb2 * gtemp4 + nzb2 * gtemp7;
                        double nz2 = nxb2 * gtemp2 + nyb2 * gtemp5 + nzb2 * gtemp8;

                        double argthetai = -(nx1 * dx + ny1 * dy + nz1 * dz);
                        double argthetaj = (nx2 * dx + ny2 * dy + nz2 * dz);

                        // cout << p1 << " " << p2 << " " << wp1 << " " << wp2 << " " << disp << " " << thetam << endl;
                        // pausel();
                        //cout << disp << endl;
                        //different conditions depending on whether there is binding or not.

                        double disp2;
                        bool cond1 = bo.boundto[wp1] == wp2 && bo.boundto[wp2] == wp1;
                        bool b1, b2;

                        b1 = bo.isbound[wp1];
                        b2 = bo.isbound[wp2];

                        if (b1 && b2 && cond1)
                        { //both bound and to each other
                            disp2 = disp;
                        }
                        else if (b1 && b2 && !cond1) //both bound and not to each other
                        {
                            disp2 = 0.5 * disp; //if both bound, make the conditions more onerous
                        }
                        else if (!b1 != !b2) //only one bound
                        {
                            disp2 = 0.7 * disp; //more onerous
                        }
                        else
                        {
                            //neither bound
                            disp2 = disp;
                        }

                        if (argthetai > ctm && argthetaj > ctm && dis < disp2)
                        {

                            double scr1 = 1 - (argthetai - cos(thetam));
                            double scr2 = 1 - (argthetaj - cos(thetam));

                            double scr3 = 2 * (dis / disp2);

                            double scr4 = -log(1E-10 + bm.calculate_score(wp1, wp2, b1 && b2 && cond1));

                            double scr = scr1 + scr2 + scr3 + scr4;

                            mdpairwd test(wp1, wp2, scr);
                            edgelist_private.push_back(test);
                        }
                    }
                }
                // else {
                //     pairschecked += fi-i;
                // }
                //pausel();
            }

            //     }
            // }

#pragma omp for schedule(static) ordered
            for (int i = 0; i < omp_get_num_threads(); i++)
            {
#pragma omp ordered
                edgelist.insert(edgelist.end(), edgelist_private.begin(), edgelist_private.end());
            }
        }
    }
    std::random_device rd; //to check whether this is too long
    std::mt19937 g(rd());

    std::shuffle(edgelist.begin(), edgelist.end(), g); // every day I'm shuffling

    // PairHistogramExtendedParallel(edgelist, boindices, boscores, tempbound);

    int np = edgelist.size();
    for(int i = 0  ; i < np ; i++) {
        int wp1 =  edgelist[i].a;
        int wp2 =  edgelist[i].b;

        double scr =  edgelist[i].scr;
        

        bool bin1 = bo.isbound[wp1];
        bool bin2 = bo.isbound[wp2];

        if(!bin1 && !bin2 && scr < 4. ) {
            bo.isbound[wp1] = true;
            bo.isbound[wp2] = true;
            bo.boundto[wp1] = wp2;
            bo.boundto[wp2] = wp1;
        }
        else{

        }
    }



   
}

#endif /* LANGEVINRFORCEUTIL_CPP */
