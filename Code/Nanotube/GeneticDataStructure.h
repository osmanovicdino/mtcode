#ifndef GENETICDATASTRUCTURE_H
#define GENETICDATASTRUCTURE_H



struct geneticcode {

    int no_types;
    vector1<int> *patch_num;
    matrix<int> *patch_pos;

    vector<vector<int> > interactions;

    double rate;

    vector1<int> *proportions;

    explicit geneticcode(vector<string> &s){ //take a vector of strings and fill in the constructor
        
        cout << "constructor called" << endl;

        int line = 0;

        no_types = stoi(s[line]);
        patch_num = new vector1<int> (no_types);
        patch_pos = new matrix<int> (no_types,6); //6 as this is the possible arrangements
        proportions = new vector1<int> (no_types);
        
        line++;

        for(int i = 0 ; i < no_types ; i++) {
            vector<int> digits;

            for (char c: s[line]) {
                digits.push_back(c-'0');
            }
            
            if(digits.size()!=7) error("error in digit size");

            (*patch_num)[i] = digits[0]; //assign patch numbers

            for(int k = 1 ; k < digits.size() ; k++) {
                 patch_pos->gpcons(i,k-1) = digits[k]; //patch positions
            }
            line++; //go up lines
            
        }

        interactions.resize(no_types * (no_types + 1) / 2);

        int ite = 0;
        for(int i = 0  ; i < no_types ; i++) {
            for(int j = i ; j < no_types ; j++ ) {
                
                int reddy = patch_num->gpcons(i)*patch_num->gpcons(j);
                interactions[ite].assign(reddy,0);
                ite++;
            }
        }


        // cout << *patch_num << endl;
        ite = 0;
        for (int i = 0; i < no_types; i++)
        {
            for (int j = i; j < no_types; j++)
            {
                
                
                
                vector<int> digits;
                for (char c : s[line])
                {
                    digits.push_back(c - '0');
                }
                // cout << s[line] << endl;
                // cout << interactions[ite].size() << " " << digits.size() << endl;
                if(interactions[ite].size() != digits.size()) error("error in interaction string length");

                for(int k = 0  ; k < digits.size() ; k++)
                interactions[ite][k] = digits[k]; //on or off
                ite++;
                line++;
            }
        }

 
        double exponent = (double)stoi(s[line]);
        exponent -=1.;
        rate = pow(10.,exponent);

        line++;

        vector<int> digitsprop;
        for (char c : s[line])
        {
            digitsprop.push_back(c - '0');
        }

        if(digitsprop.size() != no_types) error("error in proportion length");

        for(int j = 0 ; j < no_types ; j++)
            (*proportions)[j] = digitsprop[j];

    }

    geneticcode(const geneticcode &g) {
        no_types = g.no_types;

        patch_num = new vector1<int>(no_types);
        patch_pos = new matrix<int>(no_types, 6); // 6 as this is the possible arrangements
        proportions = new vector1<int>(no_types);

        for (int i = 0; i < no_types; i++)
        {
            (*patch_num)[i] =(*(g.patch_num))[i]; // assign patch numbers
            for (int k = 0; k < 6; k++)
            {
                patch_pos->gpcons(i, k - 1) = (g.patch_pos)->gpcons(i, k - 1); // patch positions
            }
            (*proportions)[i] = (*(g.proportions))[i]; // assign patch numbers
        }

        interactions.resize(no_types * (no_types+1)/2);
        // for (int i = 0; i < no_types; i++)
        // {
        //     for (int j = i; j < no_types; j++)
        //     {
        //         int reddy = patch_num->gpcons(i) * patch_num->gpcons(j);
        //         interactions[i * no_types + j].assign(reddy, 0);
        //     }
        // }
        interactions = g.interactions;
        // for (int i = 0; i < no_types; i++)
        //     {
        //         for (int j = i; j < no_types; j++)
        //         {
        //             interactions[ite] = g.interactions[ite];
        //         }
        //     }

        rate = g.rate;


    }

    template <class T>
    geneticcode& operator=(const geneticcode &g)
    {

        delete patch_num;
        delete patch_pos;
        delete proportions;

        no_types = g.no_types;

        patch_num = new vector1<int>(no_types);
        patch_pos = new matrix<int>(no_types, 6); // 6 as this is the possible arrangements
        proportions = new vector1<int>(no_types);

        for (int i = 0; i < no_types; i++)
        {
            (*patch_num)[i] = (*(g.patch_num))[i]; // assign patch numbers
            for (int k = 0; k < 6; k++)
            {
                patch_pos->gpcons(i, k - 1) = (g.patch_pos)->gpcons(i, k - 1); // patch positions
            }
            (*proportions)[i] = (*(g.proportions))[i]; // assign patch numbers
        }

        interactions.resize(no_types * (no_types + 1) / 2);
        // for (int i = 0; i < no_types; i++)
        // {
        //     for (int j = i; j < no_types; j++)
        //     {
        //         int reddy = patch_num->gpcons(i) * patch_num->gpcons(j);
        //         interactions[i * no_types + j].assign(reddy, 0);
        //     }
        // }
        interactions = g.interactions;
        
        // interactions.resize(no_types * no_types);
        // for (int i = 0; i < no_types; i++)
        // {
        //     for (int j = i; j < no_types; j++)
        //     {
        //         int reddy = patch_num->gpcons(i) * patch_num->gpcons(j);
        //         interactions[i * no_types + j].assign(reddy, 0);
        //     }
        // }

        // for (int i = 0; i < no_types; i++)
        // {
        //     for (int j = i; j < no_types; j++)
        //     {
        //         interactions[i * no_types + j] = g.interactions[i * no_types + j];
        //     }
        // }

        rate = g.rate;

        return *this;
    }

    ~geneticcode() {
        delete patch_num;
        delete patch_pos;
        delete proportions;
    }

};

#endif
