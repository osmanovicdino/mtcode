
void SingleHistogram(vector1<int> &indexes2, matrix<int> &boindices2, vector1<int> &ccs)
{

    // when the index of the group has been assigned, this function bins them
    if ((indexes2.getsize() != boindices2.getnrows()) || ((indexes2.getsize() != ccs.getsize())))
        error("initial arrays must be same size in Single Histogram");

    int sl = boindices2.getncols();

    for (int i = 0; i < indexes2.getsize(); ++i)
    {
        int wp1 = indexes2[i];
        //mtx.lock();
        //const std::lock_guard<std::mutex> lock(mtx);

        int iterator1 = ccs[wp1];
        if (iterator1 < sl)
        {

            boindices2(wp1, iterator1) = i;
            ccs[wp1]++;
            // mtx.unlock();
        }
        //mtx.unlock();
    }
}

matrix<int> SingleHistogram(vector1<int> &indexes2, vector1<int> &ccs)
{
    int n = indexes2.getsize();
    // when the index of the group has been assigned, this function bins them
    if (((n != ccs.getsize())))
        error("initial arrays must be same size in Single Histogram");


    // vector<int> arr(n);
    // for(int i = 0  ; i  < n ; i++) {
    //     arr[i] = indexes2[i];
    // }
    // sort(arr.begin(), arr.end());

    // // find the max frequency using linear traversal
    // int max_count = 1, res = arr[0], curr_count = 1;
    // for (int i = 1; i < n; i++)
    // {
    //     if (arr[i] == arr[i - 1])
    //         curr_count++;
    //     else
    //     {
    //         if (curr_count > max_count)
    //         {
    //             max_count = curr_count;
    //             res = arr[i - 1];
    //         }
    //         curr_count = 1;
    //     }
    // }

    // // If last element is most frequent
    // if (curr_count > max_count)
    // {
    //     max_count = curr_count;
    //     res = arr[n - 1];
    // }

   // cout << max_count << endl;


    //matrix<int> boindices2(n,max_count);
    //int sl = boindices2.getncols();

    for (int i = 0; i < n; ++i)
    {
        int wp1 = indexes2[i];
        //mtx.lock();
        //const std::lock_guard<std::mutex> lock(mtx);
        //int iterator1 = ccs[wp1];
        //boindices2(wp1, iterator1) = i;
        ccs[wp1]++;
            // mtx.unlock();        
        //mtx.unlock();
    }

    int sl =  maxval(ccs);
    vector1<int> counts(n);
    matrix<int> boindices2(n, sl);
    for (int i = 0; i < n; ++i)
    {
        int wp1 = indexes2[i];
        //mtx.lock();
        //const std::lock_guard<std::mutex> lock(mtx);
        int iterator1 = counts[wp1];
        boindices2(wp1, iterator1) = i;
        counts[wp1]++;
        // mtx.unlock();
        //mtx.unlock();
    }
    return boindices2;
}

template <typename T>
void PairHistogram(vector<T> &edgelist, matrix<int> &boindices, vector1<int> &tempbound)
{
    //bin an edgelist
    //Boindices must be wide enough to store the histogram, this is the responsibility of the programmer.
    for (int i = 0; i < edgelist.size(); i++)
    {
        int wp1,wp2;
        edgelist[i].get(wp1,wp2);

        boindices(wp1, tempbound[wp1]) = wp2;
        boindices(wp2, tempbound[wp2]) = wp1;
        tempbound[wp1]++;
        tempbound[wp2]++;
    }
}

void PairHistogramExtended(vector<mdpairwd> &edgelist, matrix<mdpairwd> &boindices, vector1<int> &tempbound)
{
    //bin an edgelist
    //Boindices must be wide enough to store the histogram, this is the responsibility of the programmer.
    int n =  edgelist.size();
    for (int i = 0; i < n; i++)
    {
        int wp1 = edgelist[i].a;
        int wp2 = edgelist[i].b;
        double scr = edgelist[i].scr;


        mdpairwd m1(wp2, wp2, scr);
        mdpairwd m2(wp1, wp1, scr);

        boindices(wp1, tempbound[wp1]) = m1;
        boindices(wp2, tempbound[wp2]) = m2;
        tempbound[wp1]++;
        tempbound[wp2]++;
    }
}

void PairHistogramExtended(vector<mdpairwd> &edgelist, matrix<int> &boindices, matrix<double> &boscores, vector1<int> &tempbound)
{
    //bin an edgelist
    //Boindices must be wide enough to store the histogram, this is the responsibility of the programmer.
    int n = edgelist.size();
    for (int i = 0; i < n; i++)
    {
        int wp1,wp2;
        double scr;

        edgelist[i].gets(wp1,wp2,scr);
      //  mdpairwd m1(wp2, wp2, scr);
      //  mdpairwd m2(wp1, wp1, scr);

        boindices(wp1, tempbound[wp1]) = wp2;
        boindices(wp2, tempbound[wp2]) = wp1;
        boscores(wp1, tempbound[wp1]) = scr;
        boscores(wp2, tempbound[wp2]) = scr;
        tempbound[wp1]++;
        tempbound[wp2]++;
    }
}

void DFUtil(int i, vector1<bool> &visited, matrix<int> &adj, vector1<int> &lens)
{
    visited[i] = true;

    cout << i << " ";

    for (int j = 0; j < lens[i]; j++)
    {
        int p = adj(i, j);
        if (!visited[p])
        {
            DFUtil(p, visited, adj, lens);
        }
    }
}

void DFUtilstore(int i, vector1<bool> &visited, matrix<int> &adj, vector1<int> &lens, vector1<int> &indexes, int &iter, int &len_dem)
{
    visited[i] = true;

    //cout << i << " ";
    indexes[iter] = i;
    iter++;
    len_dem++;

    for (int j = 0; j < lens[i]; j++)
    {
        int p = adj(i, j);
        if (!visited[p])
        {
            DFUtilstore(p, visited, adj, lens, indexes, iter, len_dem);
        }
    }
}

vector1<int> ConnectedComponents(matrix<int> &adj, vector1<int> &lens, vector1<int> &indexes)
{
    //find connected components of adjacency graph
    //return the vector of all the demarkers between one component and another
    //index returns

    int n = lens.getsize();

    if (indexes.getsize() != n || adj.getnrows() != n)
        error("error in sizes of input vectors in connected components");

    int iter = 0;
    int len_dem = 0;

    vector1<int> demarkus(n);

    int num_con_comp = 0;
    vector1<bool> already_accounted(n);

    for (int i = 0; i < n; i++)
    {
        if (already_accounted[i] == false)
        {
            DFUtilstore(i, already_accounted, adj, lens, indexes, iter, len_dem);

            demarkus[num_con_comp] = len_dem;

            num_con_comp++;
            len_dem = 0;

            // cout << "\n";
        }
    }

    vector1<int> nbins(num_con_comp + 1);
    nbins[0] = 0;
    int runningtotal = 0;
    for (int i = 0; i < num_con_comp; i++)
    {
        runningtotal += demarkus[i];
        nbins[i + 1] = runningtotal;
    }

    return nbins;
}

template<typename T>
vector1<int> ConnectedComponents(vector<T> &adj, vector1<int> &indexes)
{
    //find connected components of adjacency graph
    //return the vector of all the demarkers between one component and another
    //index returns
    int n = indexes.getsize();
    vector1<int> lens(n);

    int depth_mat =  n; //can change if this is too much memory
    matrix<int> boindices(n,depth_mat);
    PairHistogram(adj, boindices, lens);

    vector1<int> demarkus = ConnectedComponents(boindices,lens,indexes);

    return demarkus;

   
}

int Maximum_ConnectedComp_Size(const vector1<int>&dd) {
    int max =  dd.gpcons(1)-dd.gpcons(0);
    for(int i = 1 ; i < dd.getsize()-1 ; i++) {
        int newmax = dd.gpcons(i+1)-dd.gpcons(i);
        if( newmax > max) max = newmax;
    }
    return max;
}

inline void UpdateMin(int &lhs, int &rhs) {
    if(lhs < rhs)
        lhs = rhs;
}

/*GRAPHING ALGORITHM BASED ON THE SV ALGORITHM 
Sources: https://www.sciencedirect.com/science/article/abs/pii/0196677482900086
         https://arxiv.org/pdf/1910.05971.pdf

*/

struct mdpairsort
{
    bool operator()(const mdpair &a1, const mdpair &a2) {
        return a1.a < a2.a;
    }
};

void ConnectedComponentsParallel(matrix<int> &adj, vector1<int> &indexes) {

    //this creates a list where each index gives the component

    //Adj is the edgelist of 2 the number of edges, each edge has to be doubly connected
    //vector1<int> f(indexes);
    vector1<int> fnext(indexes);
    //if(res.size() != indexes.size) error("error in capacity of save function");
    vector1<int> ftemp(indexes.getsize());

    // #pragma omp parallel for schedule(static)
    // for(int i = 0  ; i < indexes.getNsafe() ; i++) {
    //     f[i] =  indexes[i];
    //     gf[i] = indexes[i];
    // }
 

    int nr = adj.nrows;

    for(;;) {
        bool equiv;

        #pragma omp parallel for schedule(static)
        for(int i= 0 ; i < indexes.size ; i++) {
            ftemp.data[i] = indexes.data[i];
        }
        //#pragma omp parallel 
        
           // cout << omp_get_num_threads() << endl;
        #pragma omp parallel for schedule(static)
        for(int i = 0  ; i < nr ; i++) {
            int u = adj.mat[i*2 + 0];
            int v = adj.mat[i*2 + 1];
            // int u,v;
            // sort_doublet(u1,v1,u,v);
            
            int fu = indexes.data[u];
            //cout << u << " " << v << " " << fu << " " << f.data[v] << endl;
            if (fu == indexes.data[fu] && indexes.data[v] < fu)
            {
                fnext.data[fu] = indexes.data[v];
            }
            //pausel();
        }

        //cout << "loop done" << endl;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < indexes.size; i++)
            indexes.data[i] = fnext.data[i];         

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < indexes.size; i++)
        {
            int u = i;
            int fu = indexes.data[u];

            if (fu != indexes.data[fu])
            {
                fnext.data[u] = indexes.data[fu];
            }
        }

        equiv = true;
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < indexes.size; i++)
        {
            if (indexes.data[i] != fnext.data[i])
            {
                indexes.data[i] = fnext.data[i];
            }
            if(fnext.data[i] != ftemp.data[i]) {
                equiv = false;
            }
        }
        

        //cout << "done 3" << endl;    
        //bool equiv = true;
       if(equiv) break;


    }
    

}


//for all templates Y which have member functions a and b

template <typename Y>
void ConnectedComponentsParallel(vector<Y> &adj, vector1<int> &indexes) {
    
    //WE duplicate the connectivity naturally here

    //vector1<int> f(indexes);
    vector1<int> fnext(indexes);
    //if(res.size() != indexes.size) error("error in capacity of save function");
    vector1<int> ftemp(indexes.getsize());

    // #pragma omp parallel for schedule(static)
    // for(int i = 0  ; i < indexes.getNsafe() ; i++) {
    //     f[i] =  indexes[i];
    //     gf[i] = indexes[i];
    // }
 

    int nr = adj.size();


    for(;;) {
        //cout << "iterate" << endl;
        bool equiv;

        #pragma omp parallel for schedule(static)
        for(int i= 0 ; i < indexes.size ; i++) {
            ftemp.data[i] = indexes.data[i];
        }
        //#pragma omp parallel 
        
           // cout << omp_get_num_threads() << endl;
        #pragma omp parallel for schedule(static)
        for(int i = 0  ; i < 2*nr ; i++) {
            int u,v;
            if(i%2==0) {
            u = adj[i/2].a;
            v = adj[i/2].b;
            }
            else{
            u = adj[i/2].b;
            v = adj[i/2].a;
            }
            // int u,v;
            // sort_doublet(u1,v1,u,v);
            
            int fu = indexes.data[u];
            //cout << u << " " << v << " " << fu << " " << f.data[v] << endl;
            if (fu == indexes.data[fu] && indexes.data[v] < fu)
            {
                fnext.data[fu] = indexes.data[v];
            }

            // int u1 = v;
            // int v1 = u;
            // // int u,v;
            // // sort_doublet(u1,v1,u,v);

            // int fu1 = indexes.data[u1];
            // //cout << u << " " << v << " " << fu << " " << f.data[v] << endl;
            // if (fu1 == indexes.data[fu1] && indexes.data[v1] < fu1)
            // {
            //     fnext.data[fu1] = indexes.data[v1];
            // }
            //pausel();
        }

        //cout << "loop done" << endl;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < indexes.size; i++)
            indexes.data[i] = fnext.data[i];         

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < indexes.size; i++)
        {
            int u = i;
            int fu = indexes.data[u];

            if (fu != indexes.data[fu])
            {
                fnext.data[u] = indexes.data[fu];
            }
        }

        equiv = true;
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < indexes.size; i++)
        {
            if (indexes.data[i] != fnext.data[i])
            {
                indexes.data[i] = fnext.data[i];
            }
            if(fnext.data[i] != ftemp.data[i]) {
                equiv = false;
            }
        }
        

        //cout << "done 3" << endl;    
        //bool equiv = true;
       if(equiv) break;


    }

   // pausel();
    

}

bool IndependentEdge(const vector<mdpair> &pairs, const vector1<bool> &bonds) {
    int bs = bonds.getsize();
    if(pairs.size() != bs) error("error in independent edge calculation");

    bool result = true;
    vector<int> iterators;
    //int k=0;
    iterators.reserve(bs*2);
    for(int i  = 0 ; i < bs ; i++) {
        if(bonds.gpcons(i)) {
            iterators.push_back(pairs[i].a);
            iterators.push_back(pairs[i].b);

            sort(iterators.begin(), iterators.end());
            auto it = std::unique(iterators.begin(), iterators.end());
            bool wasUnique = (it == iterators.end());
            if(!wasUnique) return false;
        }
    }

    return result;

}