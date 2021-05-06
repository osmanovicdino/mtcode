#ifndef NANOSTAR_CPP
#define NANOSTAR_CPP

Nanostar::Nanostar(matrix<double> startPositions, int N, double ll) :  placement(startPositions), bindpairs(vector<mdpair>()), bendtriples(vector<mdtriplet>()), stickerList(vector<int>()), particles(matrix<double>()){
    obj = new LangevinNVT;
    dimension = 3;
    l = ll;
    num_nanostars = N;
    total_particles = 0;
    // matrix<double> placement = startPositions;
    num_branches = 4;

    length_of_branch = 10;

    vector1<bool> is_periodic(dimension,true);
    cube bc(ll,is_periodic,dimension);

    double sigma = 1.0;

    double att_epp = 300.0; // mocified this 

    // change this line below
    // increase by a factor of 10

    WCAPotential StickerPotential(10.0,sigma,att_epp);

    WCAPotential HS(10.0, sigma, 0.0);

    FENEPotential nrr(50,1.4);

    double preferred_angle = 0.0;

    double bending_strength = 15.0;

    BendingPotential nfr2(bending_strength,preferred_angle);

    // for(int i =  0 ; i < num_nanostars; i++) {
    //     this->Passa_set_nanostar(30, 20, 4, 3, 5, "test.csv");


    // }

    int armCount = 4;
    int armLength = 3;
    setNanostar(armCount, armLength, 120, "test.csv");
    // cout << "past setnanostar" << '\n';
    // Passa_set_nanostar(start, 30, 20, 4, 3, "test.csv");
    for (int p = 0; p < num_nanostars; p++){
      int centerIdx = p*(armLength*armCount + 1);
      sortPairsTriplets(centerIdx, armCount, armLength);
      initStickerList(centerIdx, armCount, armLength);
    }

    // cout << "past sorting" <<"\n";
    potential *q1 = StickerPotential.clone();
    faa = q1;
    potential *q2 = nrr.clone();
    bindp = q2;
    potential3 *q3 = nfr2.clone();
    bendp = q3;
    potential *q4 = HS.clone();
    hs = q4;


    // matrix<double> store = create_initial_state();
    LangevinNVT b(bc);
    // cout << "after langevinnvt" <<"\n";
    double kT = 1.0;
    double dt = 0.005;
    double eta = 50.;
    //gamma = eta;

    b.setdat(particles); // <-- particles
    // call sortPairsTriplets
    //b.setinteractions(nswca);
    b.setkT(kT);
    b.setdt(dt);
    b.setgamma(eta);
    b.setm(1.0);
    matrix<double> moms(particles.getNsafe(), dimension);
    b.setmom(moms);

    *obj = b;
}

matrix<double> Nanostar::create_initial_state() {

matrix<double> store(total_particles,3);

return store;

}


matrix<double> Nanostar::create_initial_state(string s)
{

    double T;
    bool err;
    matrix<double> store = importcsv(s,T,err);
    if (store.getNsafe() !=  total_particles)
        error("size of data in filename A not correct (num of particles)");

    if (store.getncols() != dimension)
        error("size of data in filename A not correct (dimension of space)");

    if(err) error("IMPORT FILE NOT FOUND!");

    return store;
}

void Nanostar::setNanostar(int arms, int armLength, double theta, string fileName){
  int totalRealParticles = num_nanostars*(arms*armLength + 1);
  matrix<double> storeParticles(totalRealParticles, 3);
  int iter = 0;
  double pi = 2*acos(0.0);
  double incrementAngle = 2*pi / (arms - 1);
  theta = convertToRadians(theta);
  // std::ofstream myFile(fileName);
  // cout << placement << "\n";
  pausel();
  for (int j = 0; j < num_nanostars; j++){
    vector1<double> start = placement.getrowvector(j);

    //string centerCoordinate = to_string(start[0]) + ',' + to_string(start[1]) + ',' + to_string(start[2]);
    // cout << "placing center" << endl;
    storeParticles(iter,0) = start[0];
    storeParticles(iter,1) = start[1];
    storeParticles(iter,2) = start[2];
    iter++;

  //  myFile << centerCoordinate;
    for (int i = 0; i < arms - 1; i++) {
      double phi = incrementAngle*i;
      double xEnd = armLength * sin(theta)*cos(phi);
      double yEnd = armLength * sin(theta)*sin(phi);
      double zEnd = armLength * cos(theta);
      vector<double> xList = linspace(start[0], xEnd + start[0], armLength + 1);
      vector<double> yList = linspace(start[1], yEnd + start[1], armLength + 1);
      vector<double> zList = linspace(start[2], zEnd + start[2], armLength + 1);

      // placing center point
      for (int h = 1; h < armLength + 1; h++)
      {

        double x = applyPeriodicBC(xList[h], l);
        double y = applyPeriodicBC(yList[h], l);
        double z = applyPeriodicBC(zList[h], l);
      //  string coordinateToPrint = to_string(x) + ',' + to_string(y) + ',' + to_string(z) + '\n';

        //myFile << coordinateToPrint;
        // cout << j << " " << i << " " << h << endl;
        storeParticles(iter,0) = x;
        storeParticles(iter,1) = y;
        storeParticles(iter,2) = z;

        iter++;
        // check if x, y, z coords are in bounds
        // if not, apply periodic boundary condition

      }
        }
      // the last arm that comes straight up
      vector<double> zListLastArm = linspace(start[2], start[2] + armLength, armLength + 1);
      for (int d = 1; d < armLength + 1; d++)
      {
        double x = applyPeriodicBC(start[0], l);
        double y = applyPeriodicBC(start[1], l);
        double z = applyPeriodicBC(zListLastArm[d], l);
        // string coordinateToPrint = to_string(x) + ',' + to_string(y) + ',' + to_string(z) +'\n';
        // myFile << coordinateToPrint;
        cout << j << " " << 3 << " " << d << endl;
        storeParticles(iter,0) = x;
        storeParticles(iter,1) = y;
        storeParticles(iter,2) = z;
        iter++;
      }

  }
  pausel();


  // storing csv output to file

  // store = importcsv(fileName, T, err);
  // (*obj).setdat(store);
  // write to particles with storeParticles Matrix 
  particles = storeParticles;
}

void Nanostar::Passa_set_nanostar(vector1<double>start, double theta, double phi, int arms, int armLength, string fileName) {
    int totalParticles = arms * armLength;
    int totalRealParticles = totalParticles + 1; // add the nanostarCenter
    matrix <double> temp(totalRealParticles, 3);
    matrix<double> store(3,3);
    double pi = 2*acos(0.0);
 // center the box at (0, 0, 0)

    theta = convertToRadians(theta);
    phi = convertToRadians(phi);

    double incrementAngle = 2*pi / arms;

    matrix<double> I (3,3); // declaring 3x3 identity matrix
    I(0, 0) = 1;
    I(1, 1) = 1;
    I(2, 2) = 1;


    matrix<double> K (3,3); // rotating about z-axis
    K(0, 1) = -1;
    K(1, 0) = 1;

    matrix <double> K2 = K*K;

    matrix<double> R = I + sin(incrementAngle)*K + (1 - cos(incrementAngle))*K2; // Rodrigues formula

    vector1 <double> currentMaxArmCoords(3);
    currentMaxArmCoords[0] = armLength * sin(theta) * cos(phi);
    currentMaxArmCoords[1] = armLength * sin(theta) * sin(phi); // verify later
    currentMaxArmCoords[2] = armLength *cos(theta);

    // init file i/o for csv file

    std::ofstream myFile(fileName);
    for (int i = 0; i < arms; i++)
    {
      std::vector<double> xList = linspace(start[0], currentMaxArmCoords[0] + start[0], armLength + 1);
      std::vector<double> yList = linspace(start[1], currentMaxArmCoords[1] + start[1], armLength + 1);
      std::vector<double> zList = linspace(start[2], currentMaxArmCoords[2] + start[2], armLength + 1);
      for (int h = 0; h < armLength + 1; h++)
      {
        if (i != 0 && h == 0)
        {
          h = 1;
        }
        // check if x, y, z coords are in bounds
        // if not, apply periodic boundary condition
        double x = applyPeriodicBC(xList[h], l);
        double y = applyPeriodicBC(yList[h], l);
        double z = applyPeriodicBC(zList[h], l);
        string coordinateToPrint = to_string(x) + ',' + to_string(y) + ',' + to_string(z) + '\n';

        myFile << coordinateToPrint;
      }

      currentMaxArmCoords = R*currentMaxArmCoords; // problem here
      // check if generalized to non origin points
    }

    myFile.close();

    // storing csv output to file
    double T;
    bool err;
    // store = importcsv(fileName, T, err);
    // (*obj).setdat(store);
    temp = importcsv(fileName, T, err);
    particles = temp;
}

void Nanostar::sortPairsTriplets(int centerIdx, int arms, int armLength)
{
  // sorting the pairs

  // iterate through the pairs

  float nanostarCenter = centerIdx;
  float currentParticleIndex = centerIdx + 1;
  for (int i = 0; i < arms; i++)
  {
    float prevParticle = nanostarCenter;
    for (int z = 0; z < armLength; z++)
    {
        mdpair currentPair;
        currentPair.a = prevParticle;
        currentPair.b = currentParticleIndex;
        bindpairs.push_back(currentPair);
        prevParticle = currentParticleIndex;
        currentParticleIndex++;
    }


  }

  currentParticleIndex = centerIdx + 1;
  for (int i = 0; i < arms; i++)
  {
    float prevParticle = nanostarCenter;

    //this is probably redundant

    for (int z = 0; z < armLength - 1; z++){
      if (i > 0 && z == 0) {
        currentParticleIndex++;
      }
      mdtriplet currentTriplet; // different loop for out of bounds
      currentTriplet.a = prevParticle;
      currentTriplet.b = currentParticleIndex;
      currentTriplet.c = currentParticleIndex + 1;
      bendtriples.push_back(currentTriplet);
      prevParticle = currentParticleIndex;
      currentParticleIndex++;

    }
  }
}

void Nanostar::initStickerList(int centerIdx, int arms, int armLength)
{
  for (int i = 1; i < arms + 1; i++)
  {
    stickerList.push_back(centerIdx + armLength*i);
  }
}

void Nanostar::inStickerList(matrix<int> &possiblePairs, matrix <int> &specials, matrix<int> &notSpecials) // &possiblePairs, &specials &notSpecials
{
  vector<mdpair> specialVector;
  vector<mdpair> notSpecialVector;
  for (int i = 0; i < possiblePairs.getNsafe(); i++)
  {
    int firstIndex = possiblePairs(i, 0);
    int secondIndex = possiblePairs(i, 1);
    mdpair temp(firstIndex, secondIndex);
    if (count(stickerList.begin(), stickerList.end(), firstIndex) > 0 && count(stickerList.begin(), stickerList.end(), secondIndex))
    {
      specialVector.push_back(temp);
    }
    else {
      notSpecialVector.push_back(temp);
    }
  }
  matrix<int> specialsTemp(specialVector.size(), 2);
  matrix<int> notSpecialsTemp(notSpecialVector.size(), 2);

  for (int i = 0; i < specialVector.size(); i++){
    specialsTemp(i, 0) = specialVector[i].a;
    specialsTemp(i, 1) = specialVector[i].b;
  }

  for (int i = 0; i < notSpecialVector.size(); i++){
    notSpecialsTemp(i, 0) = notSpecialVector[i].a;
    notSpecialsTemp(i, 1) = notSpecialVector[i].b;
  }
  specials = specialsTemp;
  notSpecials = notSpecialsTemp;
}



void Nanostar::set_initial_state(string s)
{

    double T;
    bool err;
    matrix<double> store = importcsv(s, T, err);
    if (store.getNsafe() != total_particles)
        error("size of data in filename A not correct (num of particles)");

    if (store.getncols() != dimension)
        error("size of data in filename A not correct (dimension of space)");

    if (err)
        error("IMPORT FILE NOT FOUND!");

    (*obj).setdat(store);

}

void Nanostar::create_nanostar() {


    int initial = total_particles;

    for(int i = 0 ; i < num_branches ; i++) {
        mdpair a(initial+0,initial+length_of_branch*i+1);
        bindpairs.push_back(a);
    }

    for (int i = 0; i < num_branches; i++)
    {
        mdtriplet a(initial+0, initial+length_of_branch * i + 1, initial+length_of_branch*1+2);
        bendtriples.push_back(a);
    }

    for(int nb = 0 ; nb < num_branches ; nb++)
        for(int i = 1 ; i < length_of_branch ; i++ ) {
            mdpair a(initial + nb * length_of_branch + i, initial + nb * length_of_branch + i + 1);
            bindpairs.push_back(a);
        }

    for (int nb = 0; nb < num_branches; nb++)
        for (int i = 1; i < length_of_branch - 1; i++)
        {
            mdtriplet a(initial + nb * length_of_branch + i, initial + nb * length_of_branch + i + 1, initial + nb * length_of_branch + i + 2);
            bendtriples.push_back(a);
        }
    total_particles +=num_branches*length_of_branch+1;

}

// matrix<int> &pairs, matrix<int> &specials, matrix<int> &not_specials

// matrix<int> Nanostar::gets(vector1<double> particles){
//
//
//     //return pairs;
//     // vector<mdpair> special_pairs;
//     // vector<mdpair> not_special_pairs;
//     // for(int i = 0 ; i < pairs.getNsafe() ; i++) {
//     //     int p1 = pairs(i,0);
//     //     int p2 = pairs(i,1);
//
//     //     mdpair temp(p1,p2);
//
//     //     int which_nanostar =  floor(p1/(num_branches*length_of_branch+1));
//     //     int which_nanostar2 =  floor(p2/(num_branches*length_of_branch+1));
//
//     //     bool endp1 = ((p1 - 1) % 10 == 0);
//     //     bool endp2 = ((p2 - 1) % 10 == 0);
//
//     //     if(which_nanostar != which_nanostar2  && endp1 && endp2) {
//     //         special_pairs.push_back(temp);
//     //     }
//     //     else{
//     //         not_special_pairs.push_back(temp);
//     //     }
//
//
//     // }
//
// }



void Nanostar::run(int runtime, int every, string strbase = "")
{
    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);


    matrix<int> bp2(bindpairs.size(),2);


    matrix<int> bep2(bendtriples.size(), 3);

    //Collect all the interactions

    for(int i = 0 ; i < bindpairs.size() ; i++ ) {
        mdpair temp =  bindpairs[i];
        bp2(i, 0 ) = temp.a;
        bp2(i, 1) = temp.b;
    }

    for (int i = 0; i < bendtriples.size(); i++)
    {
        mdtriplet temp = bendtriples[i];
        bep2(i, 0) = temp.a;
        bep2(i, 1) = temp.b;
        bep2(i, 2) = temp.c;
    }




    int totalN = obj->getN();
    //int sibdiv = floor(ll/4.0);
    int ccc;
    int num = floor(l/4.);

    matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *froyo1 = obj->calculatepairs(boxes, 3.5);

    // *possible_stickers = new matrix<int>;
    matrix<int> possible_stickers;
    matrix<int> hs_pairs;
    this->inStickerList(*froyo1,possible_stickers,hs_pairs);

     //matrix<int> *hs_pairs = new matrix<int>;
    //matrix<int> hs_pairs = this->get_non_s(*froyo1);

    unsigned int i;
     for (i = 0; i < runtime; i++)
     {
       // debugging purposes
       // cout << "bind pairs" << endl;
       // cout << bp2 << endl;
       // cout << "bend triples" << endl;
       // cout << bep2 << endl;
       // cout << "possible stickers" << endl;
       // cout << possible_stickers << endl;
       // cout << "hs_pairs" << endl;
       // cout << hs_pairs << endl;

      // pausel();


         cout << i << endl;
         if (i % 25 == 0)
         {
             delete froyo1;

             // cout << "updated after: " << i << endl;
             // state = obj->getdat();
             froyo1 = obj->calculatepairs(boxes, 3.5);
             this->inStickerList(*froyo1, possible_stickers, hs_pairs);
         }

         matrix<double> F1((*obj).calculateforces(bp2, *bindp));

         // UNCOMMENT TO ADD BENDING
        matrix<double> F2((*obj).calculateforces_threebody(bep2, *bendp));

         matrix<double> F3((*obj).calculateforces(possible_stickers, *faa));

         matrix<double> F4((*obj).calculateforces(hs_pairs, *hs));

         // Debugging purposes

         // cout << "binding forces" << endl;
         // cout << F1 << endl;
         //
         // cout << "bending forces" << endl;
         // cout << F2 << endl;
         //
         // cout << "sticker forces" << endl;
         // cout << F3 << endl;
         //
         // cout << "non sticker forces" << endl;
         // cout << F4 << endl;
         // pausel();

         matrix<double> F = F1 + F2 + F3 + F4;

         matrix<double> R(totalN, dimension);
         for (int i1 = 0; i1 < totalN; i1++)
         {
             for (int j = 0; j < dimension; j++)
             {
                 R(i1, j) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
             }
         }
         if (i % every == 0)
         {

             //cout << i << endl;

             stringstream ss;

             ss << setw(number_of_digits) << setfill('0') << (i / every);

             matrix<double> pos = obj->getdat();

             string poss = "pos";
             poss = poss + strbase;

             poss += "_i=";

             string extension = ".csv";

             poss += ss.str();

             poss += extension;

             ofstream myfile;
             myfile.open(poss.c_str());

             myfile <<= pos;

             myfile.close();
         }

         (*obj).advance_mom(F, R);

         (*obj).advance_pos();
    }
}

#endif /* NANOSTAR_CPP */
