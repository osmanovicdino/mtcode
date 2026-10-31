#ifndef NANOTUBE4_CPP
#define NANOTUBE4_CPP

struct anchoringpointsharmonic
{
    double epsilon=10.;
    double anchorpointx, anchorpointy, anchorpointz;

    vector1<double> force(double x, double y, double z)
    {

        vector1<double> force(3);

        force[0] = -epsilon * (x - anchorpointx);
        force[1] = -epsilon * (y - anchorpointy);
        force[2] = -epsilon * (z - anchorpointz);
        return force;
        
    }
};

void NanotubeAssembly::run_fan_equil(int runtime, int every, geneticcode &g, string strbase)
{

    // generate boxes first
    int no_types = g.no_types;
    int ccc;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    matrix<int> boxes = obj->getgeo().generate_boxes_relationships(num, ccc);

    int Nfan = 229;

    int NN = obj->getN() + Nfan;

    // we also need to define the position of the "particle" that acts as a wall on the

    matrix<double> newdat(NN, 3);

    std::array<std::array<double, 3>, 229> fan_data =
        {{
            {-5 + ll / 2, ll / 2, 11}, {-5 + ll / 2, ll / 2, 12}, {-5 + ll / 2, ll / 2, 13}, {-5 + ll / 2, ll / 2, 14}, {-5 + ll / 2, ll / 2, 15}, {-5 + ll / 2, ll / 2, 16}, {-5 + ll / 2, ll / 2, 17}, {-5 + ll / 2, ll / 2, 18}, {-5 + ll / 2, ll / 2, 19}, {-4 + ll / 2, ll / 2, 11}, {-4 + ll / 2, ll / 2, 12}, {-4 + ll / 2, ll / 2, 13}, {-4 + ll / 2, ll / 2, 14}, {-4 + ll / 2, ll / 2, 15}, {-4 + ll / 2, ll / 2, 16}, {-4 + ll / 2, ll / 2, 17}, {-4 + ll / 2, ll / 2, 18}, {-4 + ll / 2, ll / 2, 19}, {-3 + ll / 2, ll / 2, 11}, {-3 + ll / 2, ll / 2, 12}, {-3 + ll / 2, ll / 2, 13}, {-3 + ll / 2, ll / 2, 14}, {-3 + ll / 2, ll / 2, 15}, {-3 + ll / 2, ll / 2, 16}, {-3 + ll / 2, ll / 2, 17}, {-3 + ll / 2, ll / 2, 18}, {-3 + ll / 2, ll / 2, 19}, {-2 + ll / 2, ll / 2, 11}, {-2 + ll / 2, ll / 2, 12}, {-2 + ll / 2, ll / 2, 13}, {-2 + ll / 2, ll / 2, 14}, {-2 + ll / 2, ll / 2, 15}, {-2 + ll / 2, ll / 2, 16}, {-2 + ll / 2, ll / 2, 17}, {-2 + ll / 2, ll / 2, 18}, {-2 + ll / 2, ll / 2, 19}, {-1 + ll / 2, ll / 2, 11}, {-1 + ll / 2, ll / 2, 12}, {-1 + ll / 2, ll / 2, 13}, {-1 + ll / 2, ll / 2, 14}, {-1 + ll / 2, ll / 2, 15}, {-1 + ll / 2, ll / 2, 16}, {-1 + ll / 2, ll / 2, 17}, {-1 + ll / 2, ll / 2, 18}, {-1 + ll / 2, ll / 2, 19}, {ll / 2, ll / 2, 11}, {ll / 2, ll / 2, 12}, {ll / 2, ll / 2, 13}, {ll / 2, ll / 2, 14}, {ll / 2, ll / 2, 15}, {ll / 2, ll / 2, 16}, {ll / 2, ll / 2, 17}, {ll / 2, ll / 2, 18}, {ll / 2, ll / 2, 19}, {1 + ll / 2, ll / 2, 11}, {1 + ll / 2, ll / 2, 12}, {1 + ll / 2, ll / 2, 13}, {1 + ll / 2, ll / 2, 14}, {1 + ll / 2, ll / 2, 15}, {1 + ll / 2, ll / 2, 16}, {1 + ll / 2, ll / 2, 17}, {1 + ll / 2, ll / 2, 18}, {1 + ll / 2, ll / 2, 19}, {2 + ll / 2, ll / 2, 11}, {2 + ll / 2, ll / 2, 12}, {2 + ll / 2, ll / 2, 13}, {2 + ll / 2, ll / 2, 14}, {2 + ll / 2, ll / 2, 15}, {2 + ll / 2, ll / 2, 16}, {2 + ll / 2, ll / 2, 17}, {2 + ll / 2, ll / 2, 18}, {2 + ll / 2, ll / 2, 19}, {3 + ll / 2, ll / 2, 11}, {3 + ll / 2, ll / 2, 12}, {3 + ll / 2, ll / 2, 13}, {3 + ll / 2, ll / 2, 14}, {3 + ll / 2, ll / 2, 15}, {3 + ll / 2, ll / 2, 16}, {3 + ll / 2, ll / 2, 17}, {3 + ll / 2, ll / 2, 18}, {3 + ll / 2, ll / 2, 19}, {4 + ll / 2, ll / 2, 11}, {4 + ll / 2, ll / 2, 12}, {4 + ll / 2, ll / 2, 13}, {4 + ll / 2, ll / 2, 14}, {4 + ll / 2, ll / 2, 15}, {4 + ll / 2, ll / 2, 16}, {4 + ll / 2, ll / 2, 17}, {4 + ll / 2, ll / 2, 18}, {4 + ll / 2, ll / 2, 19}, {5 + ll / 2, ll / 2, 11}, {5 + ll / 2, ll / 2, 12}, {5 + ll / 2, ll / 2, 13}, {5 + ll / 2, ll / 2, 14}, {5 + ll / 2, ll / 2, 15}, {5 + ll / 2, ll / 2, 16}, {5 + ll / 2, ll / 2, 17}, {5 + ll / 2, ll / 2, 18}, {5 + ll / 2, ll / 2, 19}, {ll / 2, -5 + ll / 2, 11}, {ll / 2, -5 + ll / 2, 12}, {ll / 2, -5 + ll / 2, 13}, {ll / 2, -5 + ll / 2, 14}, {ll / 2, -5 + ll / 2, 15}, {ll / 2, -5 + ll / 2, 16}, {ll / 2, -5 + ll / 2, 17}, {ll / 2, -5 + ll / 2, 18}, {ll / 2, -5 + ll / 2, 19}, {ll / 2, -4 + ll / 2, 11}, {ll / 2, -4 + ll / 2, 12}, {ll / 2, -4 + ll / 2, 13}, {ll / 2, -4 + ll / 2, 14}, {ll / 2, -4 + ll / 2, 15}, {ll / 2, -4 + ll / 2, 16}, {ll / 2, -4 + ll / 2, 17}, {ll / 2, -4 + ll / 2, 18}, {ll / 2, -4 + ll / 2, 19}, {ll / 2, -3 + ll / 2, 11}, {ll / 2, -3 + ll / 2, 12}, {ll / 2, -3 + ll / 2, 13}, {ll / 2, -3 + ll / 2, 14}, {ll / 2, -3 + ll / 2, 15}, {ll / 2, -3 + ll / 2, 16}, {ll / 2, -3 + ll / 2, 17}, {ll / 2, -3 + ll / 2, 18}, {ll / 2, -3 + ll / 2, 19}, {ll / 2, -2 + ll / 2, 11}, {ll / 2, -2 + ll / 2, 12}, {ll / 2, -2 + ll / 2, 13}, {ll / 2, -2 + ll / 2, 14}, {ll / 2, -2 + ll / 2, 15}, {ll / 2, -2 + ll / 2, 16}, {ll / 2, -2 + ll / 2, 17}, {ll / 2, -2 + ll / 2, 18}, {ll / 2, -2 + ll / 2, 19}, {ll / 2, -1 + ll / 2, 11}, {ll / 2, -1 + ll / 2, 12}, {ll / 2, -1 + ll / 2, 13}, {ll / 2, -1 + ll / 2, 14}, {ll / 2, -1 + ll / 2, 15}, {ll / 2, -1 + ll / 2, 16}, {ll / 2, -1 + ll / 2, 17}, {ll / 2, -1 + ll / 2, 18}, {ll / 2, -1 + ll / 2, 19}, {ll / 2, 1 + ll / 2, 11}, {ll / 2, 1 + ll / 2, 12}, {ll / 2, 1 + ll / 2, 13}, {ll / 2, 1 + ll / 2, 14}, {ll / 2, 1 + ll / 2, 15}, {ll / 2, 1 + ll / 2, 16}, {ll / 2, 1 + ll / 2, 17}, {ll / 2, 1 + ll / 2, 18}, {ll / 2, 1 + ll / 2, 19}, {ll / 2, 2 + ll / 2, 11}, {ll / 2, 2 + ll / 2, 12}, {ll / 2, 2 + ll / 2, 13}, {ll / 2, 2 + ll / 2, 14}, {ll / 2, 2 + ll / 2, 15}, {ll / 2, 2 + ll / 2, 16}, {ll / 2, 2 + ll / 2, 17}, {ll / 2, 2 + ll / 2, 18}, {ll / 2, 2 + ll / 2, 19}, {ll / 2, 3 + ll / 2, 11}, {ll / 2, 3 + ll / 2, 12}, {ll / 2, 3 + ll / 2, 13}, {ll / 2, 3 + ll / 2, 14}, {ll / 2, 3 + ll / 2, 15}, {ll / 2, 3 + ll / 2, 16}, {ll / 2, 3 + ll / 2, 17}, {ll / 2, 3 + ll / 2, 18}, {ll / 2, 3 + ll / 2, 19}, {ll / 2, 4 + ll / 2, 11}, {ll / 2, 4 + ll / 2, 12}, {ll / 2, 4 + ll / 2, 13}, {ll / 2, 4 + ll / 2, 14}, {ll / 2, 4 + ll / 2, 15}, {ll / 2, 4 + ll / 2, 16}, {ll / 2, 4 + ll / 2, 17}, {ll / 2, 4 + ll / 2, 18}, {ll / 2, 4 + ll / 2, 19}, {ll / 2, 5 + ll / 2, 11}, {ll / 2, 5 + ll / 2, 12}, {ll / 2, 5 + ll / 2, 13}, {ll / 2, 5 + ll / 2, 14}, {ll / 2, 5 + ll / 2, 15}, {ll / 2, 5 + ll / 2, 16}, {ll / 2, 5 + ll / 2, 17}, {ll / 2, 5 + ll / 2, 18}, {ll / 2, 5 + ll / 2, 19}, {-2 + ll / 2, ll / 2, 1}, {-2 + ll / 2, ll / 2, 2}, {-2 + ll / 2, ll / 2, 3}, {-2 + ll / 2, ll / 2, 4}, {-2 + ll / 2, ll / 2, 5}, {-2 + ll / 2, ll / 2, 6}, {-2 + ll / 2, ll / 2, 7}, {-2 + ll / 2, ll / 2, 8}, {-2 + ll / 2, ll / 2, 9}, {-2 + ll / 2, ll / 2, 10}, {2 + ll / 2, ll / 2, 1}, {2 + ll / 2, ll / 2, 2}, {2 + ll / 2, ll / 2, 3}, {2 + ll / 2, ll / 2, 4}, {2 + ll / 2, ll / 2, 5}, {2 + ll / 2, ll / 2, 6}, {2 + ll / 2, ll / 2, 7}, {2 + ll / 2, ll / 2, 8}, {2 + ll / 2, ll / 2, 9}, {2 + ll / 2, ll / 2, 10}, {ll / 2, -2 + ll / 2, 1}, {ll / 2, -2 + ll / 2, 2}, {ll / 2, -2 + ll / 2, 3}, {ll / 2, -2 + ll / 2, 4}, {ll / 2, -2 + ll / 2, 5}, {ll / 2, -2 + ll / 2, 6}, {ll / 2, -2 + ll / 2, 7}, {ll / 2, -2 + ll / 2, 8}, {ll / 2, -2 + ll / 2, 9}, {ll / 2, -2 + ll / 2, 10}, {ll / 2, 2 + ll / 2, 1}, {ll / 2, 2 + ll / 2, 2}, {ll / 2, 2 + ll / 2, 3}, {ll / 2, 2 + ll / 2, 4}, {ll / 2, 2 + ll / 2, 5}, {ll / 2, 2 + ll / 2, 6}, {ll / 2, 2 + ll / 2, 7}, {ll / 2, 2 + ll / 2, 8}, {ll / 2, 2 + ll / 2, 9}, {ll / 2, 2 + ll / 2, 10}
        }};

    for(int i = 0  ; i < Nfan ; i++) {
        newdat(i,0) = fan_data[i][0];
        newdat(i,1) = fan_data[i][1];
        newdat(i,2) = fan_data[i][2];
    }

    anchoringpointsharmonic a1;
    anchoringpointsharmonic a2;
    anchoringpointsharmonic a3;
    anchoringpointsharmonic a4;
    a1.anchorpointx = -2 + ll / 2.;
    a1.anchorpointy = ll / 2;
    a1.anchorpointz = 0;

    a2.anchorpointx = 2 + ll / 2.;
    a2.anchorpointy = ll / 2;
    a2.anchorpointz = 0;

    a3.anchorpointx = ll / 2.;
    a3.anchorpointy = -2 + ll / 2;
    a3.anchorpointz = 0;

    a4.anchorpointx = ll / 2.;
    a4.anchorpointy = 2 + ll / 2;
    a4.anchorpointz = 0;

    vector<int> index_anchored;
    index_anchored.push_back(189);
    index_anchored.push_back(199);
    index_anchored.push_back(209);
    index_anchored.push_back(219);



    // newdat(0, 0) = ll / 2.;
    // newdat(0, 1) = ll / 2.;
    // newdat(0, 2) = h;

    obj->initialize(newdat);

    // start with a single particle
    vector<vector<int>> indices(no_types);
    for (int i = 0; i < no_types; i++)
    {
        vector<int> a;
        indices[i] = a;
    }

    indices[0].push_back(Nfan); // add 1 particle to start

    vector<vector<int>> indices_to_add(no_types);
    for (int j = Nfan+1; j < 5000+Nfan; j++)
        indices_to_add[0].push_back(j);

    for (int i = 1; i < no_types; i++)
    {
        vector<int> a;
        indices_to_add[i] = a;
        for (int j = Nfan + i * 5000; j < Nfan+ (i + 1) * 5000; j++)
            indices_to_add[i].push_back(j);
    }

    vector1<int> counts(no_types + 1);
    counts[0] = Nfan;    // index start type 0
    counts[1] = 4999+Nfan; // index start type 1
    for (int j = 2; j < no_types + 1; j++)
    {
        counts[j] = Nfan + counts[j - 1] + 5000;
    }

    vector<int> indices_fan(Nfan);
    for(int i = 0  ; i < Nfan ; i++) {
        indices_fan[i]=i;
    }

    vector<int> indices_combine = flatten(indices);


    vector<int> index_all;
    index_all.reserve(indices_fan.size() + indices_combine.size()); // preallocate memory
    index_all.insert(index_all.end(), indices_fan.begin(), indices_fan.end());
    index_all.insert(index_all.end(), indices_combine.begin(), indices_combine.end());
    ///indices_combine.insert(indices_combine.end(),indices_fan.begin(),indices_fan.end());



    matrix<int> *pairs_onlyhs = obj->calculatepairs_parallel(boxes, index_all, 3.5); // hard sphere

    matrix<int> *pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_combine, 3.5); //hard sphere
   
    matrix<int> *pairs_fan = new matrix<int>(388,2);



    std::array<std::array<int, 2>, 388> fan_index = {{{0, 1}, {0, 9}, {1, 2}, {1, 10}, {2, 3}, {2, 11}, {3, 4}, {3, 12}, {4, 5}, {4, 13}, {5, 6}, {5, 14}, {6, 7}, {6, 15}, {7, 8}, {7, 16}, {8, 17}, {9, 10}, {9, 18}, {10, 11}, {10, 19}, {11, 12}, {11, 20}, {12, 13}, {12, 21}, {13, 14}, {13, 22}, {14, 15}, {14, 23}, {15, 16}, {15, 24}, {16, 17}, {16, 25}, {17, 26}, {18, 19}, {18, 27}, {19, 20}, {19, 28}, {20, 21}, {20, 29}, {21, 22}, {21, 30}, {22, 23}, {22, 31}, {23, 24}, {23, 32}, {24, 25}, {24, 33}, {25, 26}, {25, 34}, {26, 35}, {27, 28}, {27, 36}, {27, 198}, {28, 29}, {28, 37}, {29, 30}, {29, 38}, {30, 31}, {30, 39}, {31, 32}, {31, 40}, {32, 33}, {32, 41}, {33, 34}, {33, 42}, {34, 35}, {34, 43}, {35, 44}, {36, 37}, {36, 45}, {37, 38}, {37, 46}, {38, 39}, {38, 47}, {39, 40}, {39, 48}, {40, 41}, {40, 49}, {41, 42}, {41, 50}, {42, 43}, {42, 51}, {43, 44}, {43, 52}, {44, 53}, {45, 46}, {45, 54}, {45, 135}, {45, 144}, {46, 47}, {46, 55}, {46, 136}, {46, 145}, {47, 48}, {47, 56}, {47, 137}, {47, 146}, {48, 49}, {48, 57}, {48, 138}, {48, 147}, {49, 50}, {49, 58}, {49, 139}, {49, 148}, {50, 51}, {50, 59}, {50, 140}, {50, 149}, {51, 52}, {51, 60}, {51, 141}, {51, 150}, {52, 53}, {52, 61}, {52, 142}, {52, 151}, {53, 62}, {53, 143}, {53, 152}, {54, 55}, {54, 63}, {55, 56}, {55, 64}, {56, 57}, {56, 65}, {57, 58}, {57, 66}, {58, 59}, {58, 67}, {59, 60}, {59, 68}, {60, 61}, {60, 69}, {61, 62}, {61, 70}, {62, 71}, {63, 64}, {63, 72}, {63, 208}, {64, 65}, {64, 73}, {65, 66}, {65, 74}, {66, 67}, {66, 75}, {67, 68}, {67, 76}, {68, 69}, {68, 77}, {69, 70}, {69, 78}, {70, 71}, {70, 79}, {71, 80}, {72, 73}, {72, 81}, {73, 74}, {73, 82}, {74, 75}, {74, 83}, {75, 76}, {75, 84}, {76, 77}, {76, 85}, {77, 78}, {77, 86}, {78, 79}, {78, 87}, {79, 80}, {79, 88}, {80, 89}, {81, 82}, {81, 90}, {82, 83}, {82, 91}, {83, 84}, {83, 92}, {84, 85}, {84, 93}, {85, 86}, {85, 94}, {86, 87}, {86, 95}, {87, 88}, {87, 96}, {88, 89}, {88, 97}, {89, 98}, {90, 91}, {91, 92}, {92, 93}, {93, 94}, {94, 95}, {95, 96}, {96, 97}, {97, 98}, {99, 100}, {99, 108}, {100, 101}, {100, 109}, {101, 102}, {101, 110}, {102, 103}, {102, 111}, {103, 104}, {103, 112}, {104, 105}, {104, 113}, {105, 106}, {105, 114}, {106, 107}, {106, 115}, {107, 116}, {108, 109}, {108, 117}, {109, 110}, {109, 118}, {110, 111}, {110, 119}, {111, 112}, {111, 120}, {112, 113}, {112, 121}, {113, 114}, {113, 122}, {114, 115}, {114, 123}, {115, 116}, {115, 124}, {116, 125}, {117, 118}, {117, 126}, {118, 119}, {118, 127}, {119, 120}, {119, 128}, {120, 121}, {120, 129}, {121, 122}, {121, 130}, {122, 123}, {122, 131}, {123, 124}, {123, 132}, {124, 125}, {124, 133}, {125, 134}, {126, 127}, {126, 135}, {126, 218}, {127, 128}, {127, 136}, {128, 129}, {128, 137}, {129, 130}, {129, 138}, {130, 131}, {130, 139}, {131, 132}, {131, 140}, {132, 133}, {132, 141}, {133, 134}, {133, 142}, {134, 143}, {135, 136}, {136, 137}, {137, 138}, {138, 139}, {139, 140}, {140, 141}, {141, 142}, {142, 143}, {144, 145}, {144, 153}, {145, 146}, {145, 154}, {146, 147}, {146, 155}, {147, 148}, {147, 156}, {148, 149}, {148, 157}, {149, 150}, {149, 158}, {150, 151}, {150, 159}, {151, 152}, {151, 160}, {152, 161}, {153, 154}, {153, 162}, {153, 228}, {154, 155}, {154, 163}, {155, 156}, {155, 164}, {156, 157}, {156, 165}, {157, 158}, {157, 166}, {158, 159}, {158, 167}, {159, 160}, {159, 168}, {160, 161}, {160, 169}, {161, 170}, {162, 163}, {162, 171}, {163, 164}, {163, 172}, {164, 165}, {164, 173}, {165, 166}, {165, 174}, {166, 167}, {166, 175}, {167, 168}, {167, 176}, {168, 169}, {168, 177}, {169, 170}, {169, 178}, {170, 179}, {171, 172}, {171, 180}, {172, 173}, {172, 181}, {173, 174}, {173, 182}, {174, 175}, {174, 183}, {175, 176}, {175, 184}, {176, 177}, {176, 185}, {177, 178}, {177, 186}, {178, 179}, {178, 187}, {179, 188}, {180, 181}, {181, 182}, {182, 183}, {183, 184}, {184, 185}, {185, 186}, {186, 187}, {187, 188}, {189, 190}, {190, 191}, {191, 192}, {192, 193}, {193, 194}, {194, 195}, {195, 196}, {196, 197}, {197, 198}, {199, 200}, {200, 201}, {201, 202}, {202, 203}, {203, 204}, {204, 205}, {205, 206}, {206, 207}, {207, 208}, {209, 210}, {210, 211}, {211, 212}, {212, 213}, {213, 214}, {214, 215}, {215, 216}, {216, 217}, {217, 218}, {219, 220}, {220, 221}, {221, 222}, {222, 223}, {223, 224}, {224, 225}, {225, 226}, {226, 227}, {227, 228}}};
    for(int i = 0 ; i < 388 ; i++) {
    (*pairs_fan)(i,0) = fan_index[i][0];
    (*pairs_fan)(i,1) = fan_index[i][1];
    }

    matrix<double> F(NN, 3);
    matrix<double> Fs(NN, 3);
    matrix<double> T(NN, 3);
    matrix<double> RT(NN, 6);
    matrix<double> zeromatrix(NN, 3);


    WCAPotential wsa(3.0, 1.0, 0.0);
    FENEPotential fene(30.0,1.5);



    F = obj->calculateforces(*pairs_fan,fene);

    F += obj->calculateforces(*pairs_onlyhs, wsa); // calculate the forces due to hard sphere forces

    vector1<double> dF1 = a1.force(obj->getcoordinate(index_anchored[0], 0), obj->getcoordinate(index_anchored[0], 1), obj->getcoordinate(index_anchored[0], 2));
    vector1<double> dF2 = a2.force(obj->getcoordinate(index_anchored[1], 0), obj->getcoordinate(index_anchored[1], 1), obj->getcoordinate(index_anchored[1], 2));
    vector1<double> dF3 = a3.force(obj->getcoordinate(index_anchored[2], 0), obj->getcoordinate(index_anchored[2], 1), obj->getcoordinate(index_anchored[2], 2));
    vector1<double> dF4 = a4.force(obj->getcoordinate(index_anchored[3], 0), obj->getcoordinate(index_anchored[3], 1), obj->getcoordinate(index_anchored[3], 2));
    for(int j = 0  ; j < 3 ; j++ )
        F(index_anchored[0],j) += dF1[j];
    for (int j = 0; j < 3; j++)
        F(index_anchored[1], j) += dF2[j];
    for (int j = 0; j < 3; j++)
        F(index_anchored[2], j) += dF3[j];
    for (int j = 0; j < 3; j++)
        F(index_anchored[3], j) += dF4[j];


    obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

  

    generate_uniform_random_matrix(RT, index_all); // only generate random torques for the patchy particles

    obj->create_forces_and_torques_sphere(F, T, RT, index_all, false); // only create torques and forces for patchy particles


    // particle_adder vv;
    // vv.set_indices(indices_to_add);
    // vv.set_counts(counts);

    // vector1<int> weights;

    // weights = *(g.proportions);

    // we want to choose a
    double mu = g.rate; // can potentially set different chemical potentials

    cout << mu << endl;
    // cout << *(g.proportions) << endl;
    // pausel();
    // box_vol vol2;
    // vol2.ll = ll;
    // vol2.llx = ll;
    // vol2.lly = ll;
    // vol2.llz = obj->getcoordinate(0, 2);
    // vv.set_volume(vol2);
    // vv.set_rate(g.rate);
    // vv.set_irreducible_weights(weights);

    // we choose mu such that under an ideal gas, the total number of particles will be
    //  <N> =  V exp(\mu) where V is the size of our chamber. If we choose individual \mu_i for each
    //  particle then we have to choose individual mu. However, we need to solve it such that it still follows
    //  the ideal gas scaling for the *total* number of particles, as this is going to be the most relevant
    //  quantity in pushing the lid, that we want to keep under control.

    // so we have exp(p1 x) + exp(p2 x) up to the total number of types and we want it to be equal to exp(mu)

    // perhaps this is not what we want though, we want the reserve to be larger? This would be a bit like playing with the
    //  yes I think what I will do is choose particles non randomly according to the proportion;

    vector1<double> partialmus(no_types);

    if (no_types == 1)
    {
        // proportions do not matter
        partialmus[0] = mu;
    }
    else if (no_types == 2)
    {
        int xp1 = (*(g.proportions))[0];
        int xp2 = (*(g.proportions))[1];
        // if(xp1==xp2) { //we do not need to do anything
        //     partialmus[0] = mu;
        //     partialmus[1] = mu;
        // }
        // else{
        /*
        this code sets proportions of chemical potential, the bottom one sets proportions of abundances

    double mu2 = 0.0;
    for(int k = 0 ; k < 10 ; k++) {
        //mu2 = mu2 -  ( ( exp(xp1*mu2)+exp(xp2*mu2) -exp(mu) )/(xp1*exp(xp1*mu2)+xp2*(exp(xp2*mu2)) ) );
        mu2 = mu2 - (log((exp(xp1 * mu2) + exp(xp2 * mu2) )) - mu) / ((xp1 * exp(xp1 * mu2) + xp2 * (exp(xp2 * mu2) )) / (exp(xp1 * mu2) + exp(xp2 * mu2) ));
    }
    partialmus[0]=xp1*mu2;
    partialmus[1]=xp2*mu2;
        */
        double mu2 = log(exp(mu) / ((double)(xp1 + xp2)));
        partialmus[0] = log(double(xp1)) + mu2;
        partialmus[1] = log(double(xp2)) + mu2;
        //}
    }
    else if (no_types == 3)
    {
        int xp1 = (*(g.proportions))[0];
        int xp2 = (*(g.proportions))[1];
        int xp3 = (*(g.proportions))[2];
        // if(xp1 ==  xp2 && xp2 == xp3) {
        //     partialmus[0] = mu;
        //     partialmus[1] = mu;
        //     partialmus[2] = mu;
        // }
        // else{
        /*
        double mu2 = 0.;
        for (int k = 0; k < 10; k++)
        {
            mu2 = mu2 - (log((exp(xp1 * mu2) + exp(xp2 * mu2) + exp(xp3 * mu2))) - mu) / ((xp1 * exp(xp1 * mu2) + xp2 * (exp(xp2 * mu2) + xp3 * (exp(xp3 * mu2)))) / (exp(xp1 * mu2) + exp(xp2 * mu2) + exp(xp3 * mu2)));
            //cout << mu2 << endl;
        }
        partialmus[0] = xp1 * mu2;
        partialmus[1] = xp2 * mu2;
        partialmus[2] = xp3 * mu2;
        */
        double mu2 = log(exp(mu) / ((double)(xp1 + xp2 + xp3)));
        partialmus[0] = log(double(xp1)) + mu2;
        partialmus[1] = log(double(xp2)) + mu2;
        partialmus[2] = log(double(xp3)) + mu2;
    }
    else
    {
    }

    int MC_Move = 200;


    for (int i = 0; i < runtime; i++)
    {
        // cout << i << endl;

        obj->advancemom_halfstep(F, T, index_all);

        obj->advance_pos(index_all);

        obj->rotate(indices_combine);

        if (i % MC_Move == 0)
        {
            int type_choice = rand() % no_types;
            int ins = rand() % 2;

            if (ins == 0)
            { // insert a particle
                double rangex = ll;
                double rangey = ll;
                double rangez = ll;
                double r1 = rangex * ((double)rand() / (double)(RAND_MAX));
                double r2 = rangey * ((double)rand() / (double)(RAND_MAX));
                double r3 = rangez * ((double)rand() / (double)(RAND_MAX));

                vector1<double> myvec1(3);
                myvec1[0] = r1;
                myvec1[1] = r2;
                myvec1[2] = r3;

                // we can set the particle to be at a particular place first, then calculate the energies
                int index_which = rand() % indices_to_add[type_choice].size();
                int myindex = indices_to_add[type_choice][index_which];

                obj->set_particle(myvec1, myindex);
                double NNN = (double)indices[type_choice].size();

                double en = 0;

                for(int j = 0 ; j < indices_fan.size() ; j++)  {
                    double dis = obj->distance(myindex, indices_fan[j]);
                    if (dis < 2.)
                    {
                        en += wsa.energy(dis);
                    }
                    if (en > 100.)
                    {
                        goto energ;
                    }
                }

                for (int j = 0; j < indices_combine.size(); j++)
                {

                    double dis = obj->distance(myindex, indices_combine[j]);
                    if (dis < 2.)
                    {
                        en += obj->particle_energy(myindex, indices_combine[j], *pots);
                    }
                    if (en > 100.)
                    {
                        goto energ;
                    }
                    // if any overlaps are found the particle will be rejected and we can skip the rest
                }
            energ:

                double V = rangex * rangey * rangez;

                double mymu = partialmus[type_choice];
                // cout << "add: " <<  myindex << " " << mymu << " " << en << endl;
                double acceptance = min(1., (V / (1. + NNN)) * exp(mymu) * exp(-en));

                double r = (double)rand() / (double)RAND_MAX;

                if (r < acceptance)
                {
                    indices[type_choice].push_back(myindex); // particle is added
                    remove_at(indices_to_add[type_choice], index_which);
                }
                else
                {
                }
            }
            else
            { // remove a particle
                // choose a random particle;
                if (indices[type_choice].size() > 0)
                {
                    int index_which = (rand() % indices[type_choice].size());

                    int myindex = indices[type_choice][index_which];

                    double NNN = (double)indices[type_choice].size();

                    double rangex = ll;
                    double rangey = ll;
                    double rangez = obj->getcoordinate(0, 2);

                    double en = 0.0;
                    for (int j = 0; j < indices_fan.size(); j++)
                    {
                        if (myindex != indices_fan[j])
                        {
                            double dis = obj->distance(myindex, indices_fan[j]);
                            if (dis < 2.)
                            {
                                en += wsa.energy(dis);
                            }
                        }
                    }
                    for (int j = 0; j < indices_combine.size(); j++)
                    {
                        if (myindex != indices_combine[j])
                        {
                            double dis = obj->distance(myindex, indices_combine[j]);
                            if (dis < 2.)
                            {
                                en += wsa.energy(dis);
                                en += obj->particle_energy(myindex, indices_combine[j], *pots);
                            }
                        }
                        // if any overlaps are found the particle will be rejected and we can skip the rest
                    }
                    // that's the energy currently, we want deltaU, so the energy if it's taken away

                    double V = rangex * rangey * rangez;

                    double mymu = partialmus[type_choice];
                    // cout << "remove: " << myindex << " " << mymu << " " << en << endl;

                    double acceptance = min(1., (NNN / V) * exp(-mymu) * exp(en));

                    double r = (double)rand() / (double)RAND_MAX;

                    if (r < acceptance)
                    {
                        remove_at(indices[type_choice], index_which);

                        auto pos = std::lower_bound(indices_to_add[type_choice].begin(), indices_to_add[type_choice].end(), myindex);
                        indices_to_add[type_choice].insert(pos, myindex);
                    }
                }
            }

            indices_combine = flatten(indices); // all present particles updated
        }

        // after the particle is added, pairs need to be recomputed. This happens so long as i % MC and i % pairs are both zero

        if (i > 0 && i % 20 == 0)
        {

            // pairs = obj->calculatepairs(boxes, 3.5);

            // ADD THE PARTICLE ADDITION METHOD HERE
            // bool dd = false;
            // vector1<double> v1(3);
            // int fi;

            // cout << obj->getcoordinate(0, 2) << endl;
            // box_vol vol2;
            // vol2.ll = ll;
            // vol2.llx = ll;
            // vol2.lly = ll;
            // vol2.llz = obj->getcoordinate(0, 2) - 1;

            // // cout << "created box vol" << endl;

            // vv.set_volume(vol2);
            // vv.add_p_w(*obj, indices, dd, v1, fi);
            // cout << "added particle" << endl;
            // if (dd)
            // {

            //     // cout << "added: " << fi << endl;

            //     indices.push_back(fi);

            //     obj->set_particle(v1, fi);
            // }


            delete pairs_onlyb;
            delete pairs_onlyhs;

            index_all.clear();
            index_all.reserve(indices_fan.size() + indices_combine.size()); // preallocate memory
            index_all.insert(index_all.end(), indices_fan.begin(), indices_fan.end());
            index_all.insert(index_all.end(), indices_combine.begin(), indices_combine.end());

            pairs_onlyhs = obj->calculatepairs_parallel(boxes, index_all, 3.5);
            pairs_onlyb = obj->calculatepairs_parallel(boxes, indices_combine, 3.5);
            // cout << "done"  << endl;
        }

        F = obj -> calculateforces(*pairs_fan, fene);


        F += obj->calculateforces(*pairs_onlyhs, wsa); // calculate the forces due to hard sphere forces

        vector1<double> dF1 = a1.force(obj->getcoordinate(index_anchored[0], 0), obj->getcoordinate(index_anchored[0], 1), obj->getcoordinate(index_anchored[0], 2));
        vector1<double> dF2 = a2.force(obj->getcoordinate(index_anchored[1], 0), obj->getcoordinate(index_anchored[1], 1), obj->getcoordinate(index_anchored[1], 2));
        vector1<double> dF3 = a3.force(obj->getcoordinate(index_anchored[2], 0), obj->getcoordinate(index_anchored[2], 1), obj->getcoordinate(index_anchored[2], 2));
        vector1<double> dF4 = a4.force(obj->getcoordinate(index_anchored[3], 0), obj->getcoordinate(index_anchored[3], 1), obj->getcoordinate(index_anchored[3], 2));
        


        for (int j = 0; j < 3; j++)
            F(index_anchored[0], j) += dF1[j];
        for (int j = 0; j < 3; j++)
            F(index_anchored[1], j) += dF2[j];
        for (int j = 0; j < 3; j++)
            F(index_anchored[2], j) += dF3[j];
        for (int j = 0; j < 3; j++)
            F(index_anchored[3], j) += dF4[j];
        // cout << forcep1 << " ";

        // pausel();
        T.reset(0.0);
        // cout << *pairs_onlyb << endl;
        obj->calculate_forces_and_torques3D(*pairs_onlyb, *pots, F, T); // calculate the forces involved due to patchy

        generate_uniform_random_matrix(RT, index_all); // only generate random torques for the patchy particles

        obj->create_forces_and_torques_sphere(F, T, RT, index_all, false); // only create torques and forces for patchy particles

        obj->advancemom_halfstep(F, T, index_all);


        if (i % every == 0 && i > 0)
        {

            // cout << i << endl;

            stringstream ss;

            ss << setw(number_of_digits) << setfill('0') << (i / every);

            matrix<double> orient = obj->getorientation();
            matrix<double> pos = obj->getdat();

            string poss = "pos";
            poss = poss + strbase;
            string oris = "div";
            oris = oris + strbase;

            string elli = "or";
            elli = elli + strbase;

            poss += "_i=";
            oris += "_i=";
            // elli += "_i=";

            string extension = ".csv";

            poss += ss.str();
            oris += ss.str();
            elli += ss.str();

            poss += extension;
            oris += extension;
            elli += extension;

            string orie = "orient";
            orie += extension;

            ofstream myfile;
            myfile.open(poss.c_str());

            ofstream myfile2;
            myfile2.open(oris.c_str());

            ofstream myfile3;
            myfile3.open(elli.c_str(), std::ios_base::app);

            ofstream myfile4;
            myfile4.open(orie.c_str());

            myfile <<= pos.getrowvector(0);
            myfile << "\n";
            for (int ik = 0; ik < index_all.size(); ik++)
            {
                myfile <<= pos.getrowvector(index_all[ik]);
                myfile << "\n";
                myfile2 << index_all[ik] << endl;
            }

            myfile4 <<= obj->getorientation();

            myfile.close();
            myfile2.close();
            myfile3.close();
            myfile4.close();

            // pausel();
        }
    }
    delete pairs_onlyb;
    delete pairs_onlyhs;
    delete pairs_fan;

    // in this example, we probably want to specify the weights to be such that the total rate of each remains fixed?

    // i.e. if we have  rates as (rate*p1) and (rate*p2) and rate*p3 where p1 are the proportions we want to keep the ps proportional

    // but the complication is that each weight has to change

    // we can make all the weights point to the same reference? And then as the particle is released we update the reference?

    // we can also make it so that the length of indicies to add is the number of types * 5000 (for instance)

    // there's probably a better way of doing this

    // differential rates

    // we can just add the first element of each
}

#endif