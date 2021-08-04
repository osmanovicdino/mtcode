#ifndef ANALYZEGROWTH_CPP
#define ANALYZEGROWTH_CPP


#include <unordered_map>

matrix<int> getgrowthcurve(string dir, ComboPatch &iny, int N_Largest) 
{
    //save value of the N_largest connected components
    
    // from the directory string, get all the  binding files
    vector<string> bindfiles;
    cout << dir << endl;
    return_csv_in_dir(dir,"bind", bindfiles);

    int n = bindfiles.size(); 

    


    matrix<int> RES(n,N_Largest);

    //for each file analyze the structure
    for(int filen = 0 ; filen < n ; filen++) {
       
        int TT;
        bool vv3;
        matrix<int> bindtemp = importcsv(dir + "/" + bindfiles[filen], TT, vv3); //import my file

        int tb = bindtemp.getncols();
        BinaryBindStore bbs2;
        vector1<bool> iss(tb);
        vector1<int> ist(tb);
        for (int i = 0; i < tb; i++)
        {
            iss[i] = (bool)bindtemp(0, i);
            ist[i] = bindtemp(1, i);
        }
        bbs2.isbound = iss;
        bbs2.boundto = ist;

        vector<mdpair> edgelist; 

        
        
        //we want to get an edgelist of all the particles
        for (int b_index = 0; b_index < tb; b_index++)
        {

            bool is_bound = bool(bbs2.isbound[b_index]);

            if(is_bound) {
                //get the patches
                int wp1 = b_index;
                int wp2 = bbs2.boundto[b_index];

                //get the indices
                int p1,p2;
                iny.which_particle(wp1, wp2, p1, p2);

                mdpair test(p1,p2); //now our score is the energy
                edgelist.push_back(test);
            }
        }

        string sg = "a";
        cout << bindfiles[filen] << endl;
        cout << filen << " " << edgelist.size() << endl;
        cout << endl;
        //std::vector<mdpair> jhg(total_number_of_patches);

        
        vector<int> indexes2 = ConnectedComponentsParallel(edgelist, tb);

        unordered_map<int, size_t> count; // holds count of each encountered number
        for (int i = 0; i < tb; i++)
            count[indexes2[i]]++;

        std::vector<int> vals;
        vals.reserve(count.size());
        //connected components are now saved to indexes2;
        for(auto kv : count) {
            vals.push_back(kv.second);
        }

        #if defined(_OPENMP)
        __gnu_parallel::sort(vals.begin(), vals.end(), greater<int>());
        #else
        std::sort(vals.begin(), vals.end(), greater<int>());
        #endif

        //now we have the edgelist
        
        for(int j = 0  ; j < N_Largest ; j++) {
            if (j > vals.size()) {
                RES(filen,j) = 0;
            }
            else {
                RES(filen,j) = vals[j];
            }
            
        }
        

        //from the binary bind store, generate the 

    }

    return RES;

    // bindfiles is our vector with all the strings of the files
}

matrix<int> getgrowthcurve_distance_periodic(string dir, double l, double binding_distance, int N_Largest) {
    vector<string> posfiles;
    cout << dir << endl;
    return_csv_in_dir(dir, "pos", posfiles);

    int n = posfiles.size();

    vector1<bool> pb(3, true);
    cube geo(l, pb, 3);

    vector<mdpair> edgelist;

    matrix<int> RES(n, N_Largest);
    
    //for each file analyze the structure
    for (int filen = 0; filen < n; filen++)
    {
        cout << filen << endl;
        int TT;
        bool vv3;
        matrix<double> postemp = importcsv(dir + "/" + posfiles[filen], TT, vv3); //import my file
        int N = postemp.getNsafe();

        vector<mdpair> edgelist;

        //we want to get an edgelist of all the particles


            //bool is_bound = bool(bbs2.isbound[b_index]);

        for(int i = 0  ; i < N ; i++) {
            for(int j = i+1 ; j < N ; j++) {
                if(geo.distance_less_than(postemp,i,j,binding_distance)) 
                {
                    mdpair test(i, j);
                    edgelist.push_back(test);
                }
                
            }
        }
        
    
    vector<int> indexes2 = ConnectedComponentsParallel(edgelist, N);

    unordered_map<int, size_t> count; // holds count of each encountered number
    for (int i = 0; i < N; i++)
        count[indexes2[i]]++;

    std::vector<int> vals;
    vals.reserve(count.size());
    //connected components are now saved to indexes2;
    for (auto kv : count)
    {
        vals.push_back(kv.second);
    }

    #if defined(_OPENMP)
        __gnu_parallel::sort(vals.begin(), vals.end(), greater<int>());
    #else
        std::sort(vals.begin(), vals.end(), greater<int>());
    #endif

    //now we have the edgelist

        for (int j = 0; j < N_Largest; j++)
        {
            if (j > vals.size())
            {
                RES(filen, j) = 0;
            }
            else
            {
                RES(filen, j) = vals[j];
            }
        }
    }

    return RES;
}

#endif /* ANALYZEGROWTH_CPP */


