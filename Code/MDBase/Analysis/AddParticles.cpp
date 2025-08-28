#ifndef ADDPARTICLES_CPP
#define ADDPARTICLES_CPP

matrix<double> CreateRandomSample(matrix<double> &pos, int n, double l, int boxes)
{
    int M = pos.getnrows();
    int d = pos.getncols();

    vector1<int> dim(d);
    // for (int i = 0; i < d; i++)
    // {
    //     int ij = 1;
    //     for (int j = 0; j < i; j++)
    //     {
    //         ij *= boxes;
    //     }
    //     dim[i] = ij;
    // }
    dim[0] = SQR(boxes);
    dim[1] = boxes;
    dim[2] = 1;

    vector<int> a; //vector of occupied boxes

    double l2 = l / (double)(boxes);

    for (int i = 0; i < M; i++)
    {
        vector1<int> v1(d);
        for (int j = 0; j < d; j++)
        {
            v1[j] = floor(pos(i, j) / l2);
        }
        a.push_back(scalar(dim, v1));
    }

    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end()), a.end());

    vector<int> possible_boxes;

    for (int i = 0; i < boxes; i++)
    {
        for (int j = 0; j < boxes; j++)
        {
            for (int k = 0; k < boxes; k++)
            {
                vector1<int> v2(d);
                v2[0] = i;
                v2[1] = j;
                v2[2] = k;
                possible_boxes.push_back(scalar(dim, v2));
            }
        }
    }

    int b1 = a.size();
    int b2 = possible_boxes.size();
    std::vector<int> v(b1 + b2); // 0  0  0  0  0  0  0  0  0  0
    std::vector<int>::iterator it;

    //std::sort(a.begin(), a.begin() + b1);   //  5 10 15 20 25
    std::sort(possible_boxes.begin(), possible_boxes.end()); // 10 20 30 40 50

    it = std::set_difference(possible_boxes.begin(), possible_boxes.end(), a.begin(), a.end(), v.begin());
    //  5 15 25  0  0  0  0  0  0  0
    v.resize(it - v.begin());

    std::vector<unsigned int> indices(v.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::random_device rd; // to check whether this is too long
    std::mt19937 g(rd());

    std::shuffle(indices.begin(), indices.end(),g);

    /*
    cout << possible_boxes.size() << endl;
    cout << a.size() << endl;
    cout << v.size() << endl;
    cout << indices.size() << endl;
    pausel();

    ofstream myfile;
    myfile.open("temp.csv");
    for(int i = 0 ; i < possible_boxes.size()-1 ; i++) {
        myfile << possible_boxes[i] <<",";
    }
    myfile << possible_boxes[possible_boxes.size()-1];
    myfile << endl;
    
    for (int i = 0; i < a.size() - 1; i++)
    {
        myfile << a[i] << ",";
    }
        myfile << a[a.size() - 1];
        myfile << endl;
    
    for (int i = 0; i < v.size() - 1; i++)
    {
        myfile << v[i] << ",";
    }
    myfile << v[v.size() - 1];
    myfile << endl;
    
    for (int i = 0; i < indices.size() - 1; i++)
    {
        myfile << indices[i] << ",";
    }
        myfile << indices[indices.size() - 1];
        myfile << endl;
    */

    matrix<double> posIn(n, 3);
    for (int i = 0; i < n; i++)
    {
        int cho1 = v[indices[i]];

        int num2 = floor(cho1 / SQR(boxes));
        int num3 = floor((cho1 - num2 * SQR(boxes)) / boxes);
        int num4 = cho1 - num2 * SQR(boxes) - num3 * boxes;
        vector1<int> myvec(3);

        myvec[0] = num2;
        myvec[1] = num3;
        myvec[2] = num4;

        vector1<int> chg(3);
        for (int j = 0; j < 3; j++)
        {
            posIn(i, j) = (l2 / 2.) + l2 * myvec[j];
            chg[j] = floor(posIn(i, j) / l2);
        }

        // cout << "choicce and actual" << endl;
        // cout << dim << endl;
        // cout << cho1 << endl;
        // cout << chg << endl;
        // cout << scalar(dim,chg) << endl;
        // pausel();
    }
    return posIn;
}
#endif /* ADDPARTICLES_CPP */
