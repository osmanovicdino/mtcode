#ifndef MICROTUBULERUNS_CPP
#define MICROTUBULERUNS_CPP

void Microtubule::run(int runtime, int every)
{
    //	pausel();
    int ccc;
    int totalN = obj->getN();
    //int sibdiv = floor(ll/4.0);

    matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *froyo1 = obj->calculatepairs(boxes, pai, 3.5);
    matrix<int> *froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
    matrix<int> *froyo3 = obj->calculatepairs(boxes, pci, 3.5);
    matrix<int> *froyo4 = obj->calculatepairs(boxes, pai, pbi, 3.5);
    matrix<int> *froyo5 = obj->calculatepairs(boxes, pai, pci, 3.5);
    matrix<int> *froyo6 = obj->calculatepairs(boxes, pbi, pci, 3.5);

    matrix<int> particle_states(totalN, 1);
    for (int i = 0; i < na; i++)
    {
        particle_states(i, 0) = 0;
    }
    for (int i = na; i < na + nb; i++)
    {
        particle_states(i, 0) = 1;
    }
    for (int i = na + nb; i < na + nb + nc; i++)
    {
        particle_states(i, 0) = 2;
    }

    ofstream statefile;
    statefile.open("particle_type.csv");
    statefile <<= particle_states;
    statefile.close();

    matrix<double> state(obj->getdat()); //the state of the system

    int i;
    bool up = false;

    // vector<matrix<double> > savef1forces;
    // vector<matrix<double> > savef2forces;
    // vector<matrix<double> > savef3forces;
    // vector<matrix<double> > savef4forces;
    // vector<matrix<double> > savef5forces;
    // vector<matrix<double> > savef6forces;
    // vector<matrix<double> > savef7forces;

    // vector<matrix<double> > saveftemp1forces;
    // vector<matrix<double> > saveftemp2forces;
    // vector<matrix<double> > saveftemp3forces;

    // vector<matrix<double> > savepositions;
    // vector<vector1<int> > savebound;
    // vector<vector1<double> > saveboundalongs;

    for (i = 0; i < runtime; i++)
    {

        cout << i << endl;

        //cout << i << endl;
        //cout << (*obj).avmom() << endl;
        if (i % 25 == 0)
        {
            delete froyo1;
            delete froyo2;
            delete froyo3;
            delete froyo4;
            delete froyo5;
            delete froyo6;
            // cout << "updated after: " << i << endl;
            // state = obj->getdat();
            froyo1 = obj->calculatepairs(boxes, pai, 3.5);
            froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
            froyo3 = obj->calculatepairs(boxes, pci, 3.5);
            froyo4 = obj->calculatepairs(boxes, pai, pbi, 3.5);
            froyo5 = obj->calculatepairs(boxes, pai, pci, 3.5);
            froyo6 = obj->calculatepairs(boxes, pbi, pci, 3.5);
        }

        //cout << "pairs" << endl;
        // cout << "pairings" << endl;

        matrix<double> ftemp2(totalN, dimension), ftemp3(totalN, dimension);
        //matrix<double> angforces1(nc,dimension-1),angforces2(nc,dimension-1),angforces3(nc,dimension-1);
        // cout << "matrices initialized" << endl;

        matrix<double> F1((*obj).calculateforces(*froyo1, *faa)); //calculate the forces using the pairs as an input

        matrix<double> F2((*obj).calculateforces(*froyo2, *fbb)); //calculate the forces using the pairs as an input

        matrix<double> ftemp1((*obj).calculateforces(*froyo3, *fcc)); //calculate the forces using the pairs as an input

        matrix<double> F3((*obj).calculateforces(*froyo4, *fab)); //calculate the forces using the pairs as an input

        this->ForcesDueToPositionPL(*froyo5, ftemp2); //calculate the forces using the pairs as an input

        this->ForcesDueToPositionPL(*froyo6, ftemp3); //calculate the forces using the pairs as an input

        this->CalculateBindings(*froyo5, *froyo6);

        matrix<double> F4 = this->BindingForces();

        matrix<double> F5 = this->PositionForcesDueToAngles();

        //	cout << "active forces" << endl;
        //cout << "pos forces" << endl;
        matrix<double> F6((*obj).calculateforces(*bondpairs, *bindm));

        //cout << "bond forces" << endl;
        matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets, *bendp));

        // cout << *bendtriplets << endl;
        // pausel();

        // cout << meanmat_end(F1,na+nb) << endl;
        // cout << meanmat_end(F2,na+nb) << endl;
        // cout << meanmat_end(F3,na+nb) << endl;
        // cout << meanmat_end(F4,na+nb) << endl;
        // cout << meanmat_end(F5,na+nb) << endl;
        // cout << meanmat_end(F6,na+nb) << endl;
        // cout << meanmat_end(F7,na+nb) << endl;

        // cout << meanmat_end(ftemp1,na+nb) << endl;

        // cout << meanmat_end(ftemp2,na+nb) << endl;
        // cout << meanmat_end(ftemp3,na+nb) << endl;
        //cout << "after check matrix" << endl;

        // if(chckmatrixsize(F1,10000)) {
        // 	error("F1");
        // }
        // if(chckmatrixsize(F2,10000)) {
        // 	error("F2");
        // }
        // if(chckmatrixsize(F3,10000)) {
        // 	error("F3");
        // }
        // if(chckmatrixsize(F4,10000)) {
        // 	error("F4");
        // }

        // if(chckmatrixsize(F6,10000)) {
        // 	error("F6");
        // }
        // if(chckmatrixsize(F7,100000)) {
        // 	error("F7");
        // }
        // if(chckmatrixsize(F5,10000)) {
        // 	error("F5");
        // }
        // if(chckmatrixsize(ftemp1,10000)) {
        // 	error("ftemp1");
        // }
        // if(chckmatrixsize(ftemp2,10000)) {
        // 	error("ftemp2");
        // }
        // if(chckmatrixsize(ftemp3,10000)) {
        // 	error("ftemp3");
        // }

        matrix<double> F = ftemp1 + ftemp2 + ftemp3 + F1 + F2 + F3 + F4 + F5 + F6 + F7; //+F4+F5;
        //matrix<double> Fa = angforces1+angforces2+angforces3;
        // matrix<double> F((*obj).calculateforces(*froyo));

        // matrix<double> F2((this->PositionForcesDueToAngles_mips(100.0)));

        // matrix<double> Fa(totalN,2);
        //cout << "add forces" << endl;

        //int dd = obj->getdimension();
        matrix<double> R(totalN, dimension);
        for (int i1 = 0; i1 < totalN; i1++)
        {
            for (int j = 0; j < dimension; j++)
            {
                R(i1, j) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
            }
        }

        if (i > 0 && i % every == 0)
        {
            // for(int j = 0 ; j < na+nb ; j++) {
            // if(bound[j]>0){
            // cout << "printed" << endl;
            // cout << j << endl;
            // cout << bound[j] << endl;
            // cout << bound_along[j] << endl;
            // cout << F5(j,'r') << endl;
            // cout << F4(j,'r') << endl;
            // cout << obj->getcoordinate(j,0) << " " << obj->getcoordinate(j,1) << endl;
            // cout << obj->getcoordinate(200,0) <<  " " << obj->getcoordinate(200,1) << endl;
            // cout << obj->getcoordinate(200+L,0) <<  " " << obj->getcoordinate(200+L,1) << endl;
            // }
            // }
            // cout << F6 << endl;
            // cout << F7 << endl;

            // cout << ftemp2 << endl;
            // cout << ftemp3 << endl;

            stringstream ss2;
            // ss2 <<i/every;
            // string pairlist = "list";

            stringstream kts;
            kts << (*obj).getkT();

            // stringstream epi;
            // epi << eps;

            // stringstream epieq;
            // epieq << eqeps;

            stringstream len;
            len << l;

            string extension = "_kT=" + kts.str() + "_l=" + len.str() + ".csv";

            stringstream ss;
            ss << setfill('0') << setw(8) << (i / every);
            string filename = "x";
            filename += ss.str();
            filename += extension;

            string momname = "bind";
            momname += ss.str();
            momname += extension;

            string baname = "bind_along";
            baname += ss.str();
            baname += extension;

            ofstream myfile;
            myfile.open(filename.c_str());
            myfile <<= (*obj).getdat();
            myfile.close();

            ofstream myfile2;
            myfile2.open(momname.c_str());
            myfile2 <<= bound;
            myfile2.close();

            ofstream myfile3;
            myfile3.open(baname.c_str());
            myfile3 <<= bound_along;
            myfile3.close();
        }
        //cout << "files printed" << endl;

        (*obj).advance_mom(F, R);

        (*obj).advance_pos();
    }
}

template <typename Fun>
void Microtubule::runSV(int runtime, int every, Fun func)
{
    //	pausel();
    int ccc;
    int totalN = obj->getN();
    //int sibdiv = floor(ll/4.0);

    matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *froyo1 = obj->calculatepairs(boxes, pai, 3.5);
    matrix<int> *froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
    matrix<int> *froyo3 = obj->calculatepairs(boxes, pci, 3.5);
    matrix<int> *froyo4 = obj->calculatepairs(boxes, pai, pbi, 3.5);
    matrix<int> *froyo5 = obj->calculatepairs(boxes, pai, pci, 3.5);
    matrix<int> *froyo6 = obj->calculatepairs(boxes, pbi, pci, 3.5);

    matrix<int> particle_states(totalN, 1);
    for (int i = 0; i < na; i++)
    {
        particle_states(i, 0) = 0;
    }
    for (int i = na; i < na + nb; i++)
    {
        particle_states(i, 0) = 1;
    }
    for (int i = na + nb; i < na + nb + nc; i++)
    {
        particle_states(i, 0) = 2;
    }

    ofstream statefile;
    statefile.open("particle_type.csv");
    statefile <<= particle_states;
    statefile.close();

    matrix<double> state(obj->getdat()); //the state of the system

    int i;
    bool up = false;

    // vector<matrix<double> > savef1forces;
    // vector<matrix<double> > savef2forces;
    // vector<matrix<double> > savef3forces;
    // vector<matrix<double> > savef4forces;
    // vector<matrix<double> > savef5forces;
    // vector<matrix<double> > savef6forces;
    // vector<matrix<double> > savef7forces;

    // vector<matrix<double> > saveftemp1forces;
    // vector<matrix<double> > saveftemp2forces;
    // vector<matrix<double> > saveftemp3forces;

    // vector<matrix<double> > savepositions;
    // vector<vector1<int> > savebound;
    // vector<vector1<double> > saveboundalongs;

    for (i = 0; i < runtime; i++)
    {
        cout << i << endl;

        //cout << i << endl;
        //cout << (*obj).avmom() << endl;
        if (i % 25 == 0)
        {
            delete froyo1;
            delete froyo2;
            delete froyo3;
            delete froyo4;
            delete froyo5;
            delete froyo6;
            // cout << "updated after: " << i << endl;
            // state = obj->getdat();
            froyo1 = obj->calculatepairs(boxes, pai, 3.5);
            froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
            froyo3 = obj->calculatepairs(boxes, pci, 3.5);
            froyo4 = obj->calculatepairs(boxes, pai, pbi, 3.5);
            froyo5 = obj->calculatepairs(boxes, pai, pci, 3.5);
            froyo6 = obj->calculatepairs(boxes, pbi, pci, 3.5);
        }

        //cout << "pairs" << endl;
        // cout << "pairings" << endl;

        matrix<double> ftemp2(totalN, dimension), ftemp3(totalN, dimension);
        //matrix<double> angforces1(nc,dimension-1),angforces2(nc,dimension-1),angforces3(nc,dimension-1);
        // cout << "matrices initialized" << endl;

        matrix<double> F1((*obj).calculateforces(*froyo1, *faa)); //calculate the forces using the pairs as an input

        matrix<double> F2((*obj).calculateforces(*froyo2, *fbb)); //calculate the forces using the pairs as an input

        matrix<double> ftemp1((*obj).calculateforces(*froyo3, *fcc)); //calculate the forces using the pairs as an input

        matrix<double> F3((*obj).calculateforces(*froyo4, *fab)); //calculate the forces using the pairs as an input

        this->ForcesDueToPositionPL(*froyo5, ftemp2); //calculate the forces using the pairs as an input

        this->ForcesDueToPositionPL(*froyo6, ftemp3); //calculate the forces using the pairs as an input

        this->CalculateBindings(*froyo5, *froyo6);

        matrix<double> F4 = this->BindingForces();

        matrix<double> F5 = this->PositionForcesDueToAngles();

        //	cout << "active forces" << endl;
        //cout << "pos forces" << endl;
        matrix<double> F6((*obj).calculateforces(*bondpairs, *bindm));

        //cout << "bond forces" << endl;
        matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets, *bendp));

        //cout << "after check matrix" << endl;

        matrix<double> F = ftemp1 + ftemp2 + ftemp3 + F1 + F2 + F3 + F4 + F5 + F6 + F7; //+F4+F5;
        //matrix<double> Fa = angforces1+angforces2+angforces3;
        // matrix<double> F((*obj).calculateforces(*froyo));

        // matrix<double> F2((this->PositionForcesDueToAngles_mips(100.0)));

        // matrix<double> Fa(totalN,2);
        //cout << "add forces" << endl;

        //int dd = obj->getdimension();
        matrix<double> R(totalN, dimension);
        for (int i1 = 0; i1 < totalN; i1++)
        {
            for (int j = 0; j < dimension; j++)
            {
                R(i1, j) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
            }
        }

        if (i > 0 && i % every == 0)
        {
            // for(int j = 0 ; j < na+nb ; j++) {
            // if(bound[j]>0){
            // cout << "printed" << endl;
            // cout << j << endl;
            // cout << bound[j] << endl;
            // cout << bound_along[j] << endl;
            // cout << F5(j,'r') << endl;
            // cout << F4(j,'r') << endl;
            // cout << obj->getcoordinate(j,0) << " " << obj->getcoordinate(j,1) << endl;
            // cout << obj->getcoordinate(200,0) <<  " " << obj->getcoordinate(200,1) << endl;
            // cout << obj->getcoordinate(200+L,0) <<  " " << obj->getcoordinate(200+L,1) << endl;
            // }
            // }
            // cout << F6 << endl;
            // cout << F7 << endl;

            // cout << ftemp2 << endl;
            // cout << ftemp3 << endl;

            stringstream ss2;
            // ss2 <<i/every;
            // string pairlist = "list";

            stringstream kts;
            kts << (*obj).getkT();

            // stringstream epi;
            // epi << eps;

            // stringstream epieq;
            // epieq << eqeps;

            stringstream len;
            len << l;

            string extension = "_kT=" + kts.str() + "_l=" + len.str() + ".csv";

            stringstream ss;
            ss << setfill('0') << setw(8) << (i / every);
            string filename = "x";
            filename += ss.str();
            filename += extension;

            string momname = "bind";
            momname += ss.str();
            momname += extension;

            string baname = "bind_along";
            baname += ss.str();
            baname += extension;

            ofstream myfile;
            myfile.open(filename.c_str());
            myfile <<= (*obj).getdat();
            myfile.close();

            ofstream myfile2;
            myfile2.open(momname.c_str());
            myfile2 <<= bound;
            myfile2.close();

            ofstream myfile3;
            myfile3.open(baname.c_str());
            myfile3 <<= bound_along;
            myfile3.close();
        }
        //cout << "files printed" << endl;

        (*obj).advance_mom_spatial_dependence(F, R, func);

        //cout << "momenta advanced" << endl;
        // if(trap(bound,1.0)>0) {
        // 	cout << "after momenta update" << endl;
        // 	cout << bound << endl;
        // 	cout << bound_along << endl;

        // 	pausel();

        // }

        (*obj).advance_pos();
        //cout << "positions advanced" << endl;
        // if(trap(bound,1.0)>0) {
        // 	cout << "after position update" << endl;
        // 	cout << bound << endl;
        // 	cout << bound_along << endl;

        // 	pausel();

        // }
        //this->UpdateBoundAlong();

        //	this->UpdateAngles(Fa,R2);
        //cout << "angles updated" << endl;
    }
}

template <typename Fun>
void Microtubule::runPV(int runtime, int every, Fun func)
{
    //	pausel();
    int ccc;
    int totalN = obj->getN();
    //int sibdiv = floor(ll/4.0);

    matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num, ccc);

    matrix<int> *froyo1 = obj->calculatepairs(boxes, pai, 3.5);
    matrix<int> *froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
    matrix<int> *froyo3 = obj->calculatepairs(boxes, pci, 3.5);
    matrix<int> *froyo4 = obj->calculatepairs(boxes, pai, pbi, 3.5);
    matrix<int> *froyo5 = obj->calculatepairs(boxes, pai, pci, 3.5);
    matrix<int> *froyo6 = obj->calculatepairs(boxes, pbi, pci, 3.5);

    matrix<int> particle_states(totalN, 1);
    for (int i = 0; i < na; i++)
    {
        particle_states(i, 0) = 0;
    }
    for (int i = na; i < na + nb; i++)
    {
        particle_states(i, 0) = 1;
    }
    for (int i = na + nb; i < na + nb + nc; i++)
    {
        particle_states(i, 0) = 2;
    }

    ofstream statefile;
    statefile.open("particle_type.csv");
    statefile <<= particle_states;
    statefile.close();

    matrix<double> state(obj->getdat()); //the state of the system

    int i;
    bool up = false;

    // vector<matrix<double> > savef1forces;
    // vector<matrix<double> > savef2forces;
    // vector<matrix<double> > savef3forces;
    // vector<matrix<double> > savef4forces;
    // vector<matrix<double> > savef5forces;
    // vector<matrix<double> > savef6forces;
    // vector<matrix<double> > savef7forces;

    // vector<matrix<double> > saveftemp1forces;
    // vector<matrix<double> > saveftemp2forces;
    // vector<matrix<double> > saveftemp3forces;

    // vector<matrix<double> > savepositions;
    // vector<vector1<int> > savebound;
    // vector<vector1<double> > saveboundalongs;

    for (i = 0; i < runtime; i++)
    {
        //cout << i << endl;

        cout << i << endl;
        //cout << (*obj).avmom() << endl;
        if (i % 25 == 0)
        {
            delete froyo1;
            delete froyo2;
            delete froyo3;
            delete froyo4;
            delete froyo5;
            delete froyo6;
            // cout << "updated after: " << i << endl;
            // state = obj->getdat();
            froyo1 = obj->calculatepairs(boxes, pai, 3.5);
            froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
            froyo3 = obj->calculatepairs(boxes, pci, 3.5);
            froyo4 = obj->calculatepairs(boxes, pai, pbi, 3.5);
            froyo5 = obj->calculatepairs(boxes, pai, pci, 3.5);
            froyo6 = obj->calculatepairs(boxes, pbi, pci, 3.5);
        }

        //cout << "pairs" << endl;
        // cout << "pairings" << endl;

        matrix<double> ftemp2(totalN, dimension), ftemp3(totalN, dimension);
        //matrix<double> angforces1(nc,dimension-1),angforces2(nc,dimension-1),angforces3(nc,dimension-1);
        // cout << "matrices initialized" << endl;

        matrix<double> F1((*obj).calculateforces(*froyo1, *faa)); //calculate the forces using the pairs as an input

        matrix<double> F2((*obj).calculateforces(*froyo2, *fbb)); //calculate the forces using the pairs as an input

        matrix<double> ftemp1((*obj).calculateforces(*froyo3, *fcc)); //calculate the forces using the pairs as an input

        matrix<double> F3((*obj).calculateforces(*froyo4, *fab)); //calculate the forces using the pairs as an input

        this->ForcesDueToPositionPL(*froyo5, ftemp2); //calculate the forces using the pairs as an input

        this->ForcesDueToPositionPL(*froyo6, ftemp3); //calculate the forces using the pairs as an input

        this->CalculateBindings(*froyo5, *froyo6);

        matrix<double> F4 = this->BindingForces();

        matrix<double> F5 = this->PositionForcesDueToAngles();

        //	cout << "active forces" << endl;
        //cout << "pos forces" << endl;
        matrix<double> F6((*obj).calculateforces(*bondpairs, *bindm));

        //cout << "bond forces" << endl;
        matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets, *bendp));

        //cout << "after check matrix" << endl;

        matrix<double> F = ftemp1 + ftemp2 + ftemp3 + F1 + F2 + F3 + F4 + F5 + F6 + F7; //+F4+F5;
        //matrix<double> Fa = angforces1+angforces2+angforces3;
        // matrix<double> F((*obj).calculateforces(*froyo));

        // matrix<double> F2((this->PositionForcesDueToAngles_mips(100.0)));

        // matrix<double> Fa(totalN,2);
        //cout << "add forces" << endl;

        //int dd = obj->getdimension();
        matrix<double> R(totalN, dimension);
        for (int i1 = 0; i1 < totalN; i1++)
        {
            for (int j = 0; j < dimension; j++)
            {
                R(i1, j) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
            }
        }

        if (i > 0 && i % every == 0)
        {
            // for(int j = 0 ; j < na+nb ; j++) {
            // if(bound[j]>0){
            // cout << "printed" << endl;
            // cout << j << endl;
            // cout << bound[j] << endl;
            // cout << bound_along[j] << endl;
            // cout << F5(j,'r') << endl;
            // cout << F4(j,'r') << endl;
            // cout << obj->getcoordinate(j,0) << " " << obj->getcoordinate(j,1) << endl;
            // cout << obj->getcoordinate(200,0) <<  " " << obj->getcoordinate(200,1) << endl;
            // cout << obj->getcoordinate(200+L,0) <<  " " << obj->getcoordinate(200+L,1) << endl;
            // }
            // }
            // cout << F6 << endl;
            // cout << F7 << endl;

            // cout << ftemp2 << endl;
            // cout << ftemp3 << endl;

            stringstream ss2;
            // ss2 <<i/every;
            // string pairlist = "list";

            stringstream kts;
            kts << (*obj).getkT();

            // stringstream epi;
            // epi << eps;

            // stringstream epieq;
            // epieq << eqeps;

            stringstream len;
            len << l;

            string extension = "_kT=" + kts.str() + "_l=" + len.str() + ".csv";

            stringstream ss;
            ss << setfill('0') << setw(8) << (i / every);
            string filename = "x";
            filename += ss.str();
            filename += extension;

            string momname = "bind";
            momname += ss.str();
            momname += extension;

            string baname = "bind_along";
            baname += ss.str();
            baname += extension;

            ofstream myfile;
            myfile.open(filename.c_str());
            myfile <<= (*obj).getdat();
            myfile.close();

            ofstream myfile2;
            myfile2.open(momname.c_str());
            myfile2 <<= bound;
            myfile2.close();

            ofstream myfile3;
            myfile3.open(baname.c_str());
            myfile3 <<= bound_along;
            myfile3.close();
        }
        //cout << "files printed" << endl;

        (*obj).advance_mom_particle_dependence(F, R, func);

        //cout << "momenta advanced" << endl;
        // if(trap(bound,1.0)>0) {
        // 	cout << "after momenta update" << endl;
        // 	cout << bound << endl;
        // 	cout << bound_along << endl;

        // 	pausel();

        // }

        (*obj).advance_pos();
        //cout << "positions advanced" << endl;
        // if(trap(bound,1.0)>0) {
        // 	cout << "after position update" << endl;
        // 	cout << bound << endl;
        // 	cout << bound_along << endl;

        // 	pausel();

        // }
        //this->UpdateBoundAlong();

        //	this->UpdateAngles(Fa,R2);
        //cout << "angles updated" << endl;
    }
}

template <typename Fun>
void Microtubule::runMTONLY(int runtime, int every, Fun func)
{
    //	pausel();

    int ccc;
    int totalN = obj->getN();
    //int sibdiv = floor(ll/4.0);

    matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num, ccc);

    cout << 1 << endl;
    // matrix<int> *froyo1 = obj->calculatepairs(boxes, pai, 3.5);
    // matrix<int> *froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
    matrix<int> *froyo3 = obj->calculatepairs(boxes, pci, 3.5);

    cout << 2 << endl;
    // matrix<int> *froyo4 = obj->calculatepairs(boxes, pai, pbi, 3.5);
    // matrix<int> *froyo5 = obj->calculatepairs(boxes, pai, pci, 3.5);
    // matrix<int> *froyo6 = obj->calculatepairs(boxes, pbi, pci, 3.5);

    matrix<int> *bindpairs2 = new matrix<int>(1, 1); //additional functions for the forces
    matrix<int> *bendpairs2 = new matrix<int>(1, 1); //additional functions for the forces

    if (number_of_microtubules == 2)
    {
        matrix<int> temp1(L, 2);

        for (int i = 0; i < L; i++)
        {
            temp1(i, 0) = na + nb + i;
            temp1(i, 1) = na + nb + L + i;
        }
        matrix<int> temp2;
        *bindpairs2 = temp1;
        *bendpairs2 = temp2;
    }
    else if (number_of_microtubules > 2)
    {

        matrix<int> temp1(L * (number_of_microtubules - 1), 2);
        matrix<int> temp2(L * (number_of_microtubules - 2), 3);

        int iter = 0;
        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < number_of_microtubules - 1; j++)
            {
                //cout << iter << endl;
                temp1(iter, 0) = na + nb + j * L + i;
                temp1(iter, 1) = na + nb + (j + 1) * L + i;
                iter++;
            }
        }
        int iter2 = 0;
        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < number_of_microtubules - 2; j++)
            {
                temp2(iter2, 0) = na + nb + j * L + i;
                temp2(iter2, 1) = na + nb + (j + 1) * L + i;
                temp2(iter2, 2) = na + nb + (j + 2) * L + i;
                iter2++;
            }
        }

        *bindpairs2 = temp1;
        *bendpairs2 = temp2;
    }
    else
    {
    }

    cout << "all set up" << endl;

    int i;
    bool up = false;

    // vector<matrix<double> > savef1forces;
    // vector<matrix<double> > savef2forces;
    // vector<matrix<double> > savef3forces;
    // vector<matrix<double> > savef4forces;
    // vector<matrix<double> > savef5forces;
    // vector<matrix<double> > savef6forces;
    // vector<matrix<double> > savef7forces;

    // vector<matrix<double> > saveftemp1forces;
    // vector<matrix<double> > saveftemp2forces;
    // vector<matrix<double> > saveftemp3forces;

    // vector<matrix<double> > savepositions;
    // vector<vector1<int> > savebound;
    // vector<vector1<double> > saveboundalongs;

    for (i = 0; i < runtime; i++)
    {

        cout << i << endl;

        //cout << i << endl;
        //cout << (*obj).avmom() << endl;
        // if (i % 25 == 0)
        // {

        // 	delete froyo3;

        // 	// cout << "updated after: " << i << endl;
        // 	// state = obj->getdat();
        // 	froyo1 = obj->calculatepairs(boxes, pai, 3.5);
        // 	froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
        // 	froyo3 = obj->calculatepairs(boxes, pci, 3.5);
        // 	froyo4 = obj->calculatepairs(boxes, pai, pbi, 3.5);
        // 	froyo5 = obj->calculatepairs(boxes, pai, pci, 3.5);
        // 	froyo6 = obj->calculatepairs(boxes, pbi, pci, 3.5);
        // }

        //cout << "pairs" << endl;
        // cout << "pairings" << endl;

        // matrix<double> ftemp2(totalN, dimension), ftemp3(totalN, dimension);
        // //matrix<double> angforces1(nc,dimension-1),angforces2(nc,dimension-1),angforces3(nc,dimension-1);
        // // cout << "matrices initialized" << endl;

        // matrix<double> F1((*obj).calculateforces(*froyo1, *faa)); //calculate the forces using the pairs as an input

        // matrix<double> F2((*obj).calculateforces(*froyo2, *fbb)); //calculate the forces using the pairs as an input

        matrix<double> ftemp1((*obj).calculateforces(*froyo3, *fcc)); //calculate the forces using the pairs as an input

        // matrix<double> F3((*obj).calculateforces(*froyo4, *fab)); //calculate the forces using the pairs as an input

        // this->ForcesDueToPositionPL(*froyo5, ftemp2); //calculate the forces using the pairs as an input

        // this->ForcesDueToPositionPL(*froyo6, ftemp3); //calculate the forces using the pairs as an input

        // this->CalculateBindings(*froyo5, *froyo6);

        // matrix<double> F4 = this->BindingForces();

        //matrix<double> F5 = this->PositionForcesDueToAngles();

        matrix<double> F5 = this->constantMTforce();

        //	cout << "active forces" << endl;
        //cout << "pos forces" << endl;
        matrix<double> F6((*obj).calculateforces(*bondpairs, *bindm));

        //cout << "bond forces" << endl;
        matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets, *bendp));

        matrix<double> F8(totalN, dimension);
        matrix<double> F9(totalN, dimension);

        if (number_of_microtubules == 2)
        {
            F8 = ((*obj).calculateforces(*bindpairs2, *bindm));
        }
        else if (number_of_microtubules > 2)
        {
            F8 = ((*obj).calculateforces(*bindpairs2, *bindm));
            //cout << "bond forces" << endl;
            F9 = ((*obj).calculateforces_threebody(*bendpairs2, *bendp));
        }

        // cout << *bendtriplets << endl;
        // pausel();

        // cout << meanmat_end(F1,na+nb) << endl;
        // cout << meanmat_end(F2,na+nb) << endl;
        // cout << meanmat_end(F3,na+nb) << endl;
        // cout << meanmat_end(F4,na+nb) << endl;
        // cout << meanmat_end(F5,na+nb) << endl;
        // cout << meanmat_end(F6,na+nb) << endl;
        // cout << meanmat_end(F7,na+nb) << endl;

        // cout << meanmat_end(ftemp1,na+nb) << endl;

        // cout << meanmat_end(ftemp2,na+nb) << endl;
        // cout << meanmat_end(ftemp3,na+nb) << endl;
        //cout << "after check matrix" << endl;

        // if(chckmatrixsize(F1,10000)) {
        // 	error("F1");
        // }
        // if(chckmatrixsize(F2,10000)) {
        // 	error("F2");
        // }
        // if(chckmatrixsize(F3,10000)) {
        // 	error("F3");
        // }
        // if(chckmatrixsize(F4,10000)) {
        // 	error("F4");
        // }

        // if(chckmatrixsize(F6,10000)) {
        // 	error("F6");
        // }
        // if(chckmatrixsize(F7,100000)) {
        // 	error("F7");
        // }
        // if(chckmatrixsize(F5,10000)) {
        // 	error("F5");
        // }
        // if(chckmatrixsize(ftemp1,10000)) {
        // 	error("ftemp1");
        // }
        // if(chckmatrixsize(ftemp2,10000)) {
        // 	error("ftemp2");
        // }
        // if(chckmatrixsize(ftemp3,10000)) {
        // 	error("ftemp3");
        // }

        matrix<double> F = ftemp1 + F5 + F6 + F7 + F8 + F9; //+F4+F5;

        // double maxsize = 100000.;
        // if(chckmatrixsize(ftemp1,maxsize)) {
        // 	error("ftemp1");
        // }
        // if(chckmatrixsize(F5,maxsize)) {
        // 	error("F5");
        // }
        // if (chckmatrixsize(F6, maxsize))
        // {
        // 	error("F6");
        // }
        // if (chckmatrixsize(F7, maxsize))
        // {
        // 	error("F7");
        // }
        // if (chckmatrixsize(F8, maxsize))
        // {
        // 	error("F8");
        // }
        // if (chckmatrixsize(F9, maxsize))
        // {
        // 	error("F9");
        // }
        //matrix<double> Fa = angforces1+angforces2+angforces3;
        // matrix<double> F((*obj).calculateforces(*froyo));

        // matrix<double> F2((this->PositionForcesDueToAngles_mips(100.0)));

        // matrix<double> Fa(totalN,2);
        //cout << "add forces" << endl;

        //int dd = obj->getdimension();
        matrix<double> R(totalN, dimension);
        for (int i1 = 0; i1 < totalN; i1++)
        {
            for (int j = 0; j < dimension; j++)
            {
                R(i1, j) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
            }
        }

        if (i > 0 && i % every == 0)
        {
            // for(int j = 0 ; j < na+nb ; j++) {
            // if(bound[j]>0){
            // cout << "printed" << endl;
            // cout << j << endl;
            // cout << bound[j] << endl;
            // cout << bound_along[j] << endl;
            // cout << F5(j,'r') << endl;
            // cout << F4(j,'r') << endl;
            // cout << obj->getcoordinate(j,0) << " " << obj->getcoordinate(j,1) << endl;
            // cout << obj->getcoordinate(200,0) <<  " " << obj->getcoordinate(200,1) << endl;
            // cout << obj->getcoordinate(200+L,0) <<  " " << obj->getcoordinate(200+L,1) << endl;
            // }
            // }
            // cout << F6 << endl;
            // cout << F7 << endl;

            // cout << ftemp2 << endl;
            // cout << ftemp3 << endl;

            stringstream ss2;
            // ss2 <<i/every;
            // string pairlist = "list";

            stringstream kts;
            kts << (*obj).getkT();

            // stringstream epi;
            // epi << eps;

            // stringstream epieq;
            // epieq << eqeps;

            stringstream len;
            len << l;

            string extension = "_kT=" + kts.str() + "_l=" + len.str() + ".csv";

            stringstream ss;
            ss << setfill('0') << setw(8) << (i / every);
            string filename = "x";
            filename += ss.str();
            filename += extension;

            string momname = "bind";
            momname += ss.str();
            momname += extension;

            string baname = "bind_along";
            baname += ss.str();
            baname += extension;

            ofstream myfile;
            myfile.open(filename.c_str());
            myfile <<= (*obj).getdat();
            myfile.close();

            ofstream myfile2;
            myfile2.open(momname.c_str());
            myfile2 <<= bound;
            myfile2.close();

            ofstream myfile3;
            myfile3.open(baname.c_str());
            myfile3 <<= bound_along;
            myfile3.close();
        }
        //cout << "files printed" << endl;

        (*obj).advance_mom_spatial_dependence(F, R, func);

        //cout << "momenta advanced" << endl;
        // if(trap(bound,1.0)>0) {
        // 	cout << "after momenta update" << endl;
        // 	cout << bound << endl;
        // 	cout << bound_along << endl;

        // 	pausel();

        // }

        (*obj).advance_pos();
    }
}

template <typename Fun>
void Microtubule::runMTONLY_initialstate(int runtime, int every, Fun func, double up_to_length, vector1<int> ori1, vector1<int> ori2)
{
    //ori 1 and ori2 are the directions of the force
    if(ori1.getsize()!=ori2.getsize()) error("list of particles along which average force is defined needs to be the same size");
    //	pausel();
    //in this function, create bonds arising from the initial state
    if (up_to_length > 3.5)
        error("defining too many bonds, speak to Dino");
    int ccc;
    int totalN = obj->getN();
    //int sibdiv = floor(ll/4.0);

    matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num, ccc);

    cout << 1 << endl;
    // matrix<int> *froyo1 = obj->calculatepairs(boxes, pai, 3.5);
    // matrix<int> *froyo2 = obj->calculatepairs(boxes, pbi, 3.5);
    matrix<int> *froyo3 = obj->calculatepairs(boxes, pci, 3.5);

    double epsilon = 100.0; //strength of the harmonic bond; //MAKE THIS STRONGER TO STIFFEN THE SYSTEM

    //for the intial dat
    int totop = 0;
    matrix<vector1<double>> params(totalN, totalN); //store all the parameters
    vector<mdpair> savepairs;

    for (int i = 0; i < froyo3->getNsafe(); i++)
    {
        int p1 = (*froyo3)(i, 0);
        int p2 = (*froyo3)(i, 1);

        double dis = obj->distance(p1, p2);

        if (dis < up_to_length)
        {
            vector1<double> temp(2);
            temp[0] = epsilon;
            temp[1] = dis;
            params(p1, p2) = temp;
            params(p2, p1) = temp;

            mdpair kj(p1, p2);
            savepairs.push_back(kj);
            totop++;
        }
    }

    matrix<int> temp2(totop, 2);

    for (int i = 0; i < totop; i++)
    {
        temp2(i, 0) = savepairs[i].a;
        temp2(i, 1) = savepairs[i].b;
    }

    matrix<int> *bindpairs = new matrix<int>(1, 1);
    *bindpairs = temp2;

    cout << 2 << endl;



    cout << "all set up" << endl;

    int i;
    bool up = false;

    for (i = 0; i < runtime; i++)
    {

        cout << i << endl;

        matrix<double> ftemp1((*obj).calculateforces(*froyo3, *fcc)); //calculate the forces using the pairs as an input

        matrix<double> F5 = this->constantMTforce(ori1,ori2);

        //	cout << "active forces" << endl;
        //cout << "pos forces" << endl;
        matrix<double> F6((*obj).calculateforces_sp(*bindpairs, *bindp, params));

        matrix<double> F = ftemp1 + F5 + F6; //+F4+F5;

        matrix<double> R(totalN, dimension);
        for (int i1 = 0; i1 < totalN; i1++)
        {
            for (int j = 0; j < dimension; j++)
            {
                R(i1, j) = (3.464101615 * ((double)rand() / (RAND_MAX)) - 1.732050808);
            }
        }

        if (i > 0 && i % every == 0)
        {

            stringstream ss2;
            // ss2 <<i/every;
            // string pairlist = "list";

            stringstream kts;
            kts << (*obj).getkT();

            // stringstream epi;
            // epi << eps;

            // stringstream epieq;
            // epieq << eqeps;

            stringstream len;
            len << l;

            string extension = "_kT=" + kts.str() + "_l=" + len.str() + ".csv";

            stringstream ss;
            ss << setfill('0') << setw(8) << (i / every);
            string filename = "x";
            filename += ss.str();
            filename += extension;

            string momname = "bind";
            momname += ss.str();
            momname += extension;

            string baname = "bind_along";
            baname += ss.str();
            baname += extension;

            ofstream myfile;
            myfile.open(filename.c_str());
            myfile <<= (*obj).getdat();
            myfile.close();

            ofstream myfile2;
            myfile2.open(momname.c_str());
            myfile2 <<= bound;
            myfile2.close();

            ofstream myfile3;
            myfile3.open(baname.c_str());
            myfile3 <<= bound_along;
            myfile3.close();
        }
        //cout << "files printed" << endl;

        (*obj).advance_mom_spatial_dependence(F, R, func);

        (*obj).advance_pos();
    }
}

#endif /* MICROTUBULERUNS_CPP */
