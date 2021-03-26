#ifndef GENERALPATCH_CPP
#define GENERALPATCH_CPP

// we define here a general patch model such that we may quickly be able to create new simulations of whichever
// geometry we are interested in

GeneralPatch::GeneralPatch(vector1<int> no_patches_per_typee, vector1<int> num_per_typee, matrix<double> &paramss, matrix<double> &orientt, bool flatt  = false) : ComboPatch(paramss.getnrows()), params(paramss), orient(orientt), no_patches_per_type(no_patches_per_typee), num_per_type(num_per_typee), total_patches_per_type(vector1<int>(no_patches_per_typee.getsize())), pot_starters(matrix<int>(num_per_typee.getsize(), num_per_typee.getsize()))
{
    no_types = no_patches_per_type.getsize();
    flat = flatt;


    //if(params.getnrows() != no_types*(no_types)/2 ) error("param matrix not of the correct size");
    total_patches_per_type[0] = num_per_type[0]*no_patches_per_type[0];
    for(int i = 1 ; i < no_types ; i++ ) {
        total_patches_per_type[i] = total_patches_per_type[i-1] + (num_per_type[i]-num_per_type[i-1])*no_patches_per_type[i];
    }

    int total_no_of_orients =  0;
    vector1<int> start_and_end_points(no_types+1);
    start_and_end_points[0] = 0;
    for(int i = 0  ; i < no_patches_per_type.getsize() ; i++) {
        total_no_of_orients += no_patches_per_type[i];
        start_and_end_points[i+1] =  total_no_of_orients;
    }


    if(orient.getnrows() != total_no_of_orients) error("incorrect no. of orientations");


    int ht =  no_types*(no_types+1)/2;
    i1 = new int *[ht]; // allocate an array of 10 int pointers — these are our rows

    
    int hj = 0;
    int nn = 0;
    for(int i = 0 ; i < no_types ; i++ ) {
        for(int j = i ; j< no_types; j++) {
            pot_starters(i,j) = hj;
            i1[nn] = new int [1+no_patches_per_type[i]*no_patches_per_type[j]];


            i1[nn][0] = no_patches_per_type[i]*no_patches_per_type[j];
           
            

            for (int k = 0; k < no_patches_per_type[i] * no_patches_per_type[j] ; k++) {
                i1[nn][k+1] =  hj;
                hj++;
            }
            nn++;
        }
    }

    p = &i1[0];



    safe = false;

    int iter  = 0 ;

    


    for(int i = 0 ; i < no_types ; i++) {
        for(int j =i ; j < no_types ; j++) {
            // for(int i1  = 0 ; i1 < no_patches_per_type[i] ; i1++ ) {
            //     for(int j1 = 0 ; j1 < no_pathes_per_type[j] ; j1++) {

                    int starti = start_and_end_points[i];
                    int endi = start_and_end_points[i+1];
                    int startj = start_and_end_points[j];
                    int endj = start_and_end_points[j+1];

                    for(int i2 =  starti ; i2 < endi ; i2++) {
                        for(int j1 = startj  ; j1 < endj ; j1++) {
                            potentialtheta3D *pot1;
                            if(flat)
                                pot1 = new KernFrenkelOnePatchFlatBottom(orient(i2, 0), orient(i2, 1), orient(i2, 2), orient(j1, 0), orient(j1, 1), orient(j1, 2), params(iter, 0), params(iter, 1), params(iter, 2), 0.5);
                            else 
                                pot1 = new mypot(orient(i2, 0), orient(i2, 1), orient(i2, 2), orient(j1, 0), orient(j1, 1), orient(j1, 2), params(iter, 0), params(iter, 1), params(iter, 2), 0.2);


                            potential_bundle[iter] = pot1->clone();
                            delete pot1;
                            iter++;     
                }
            }
        }
    }

    



}

GeneralPatch::GeneralPatch(const GeneralPatch &a) : ComboPatch((a.params).getnrows()), params(a.params), orient(a.orient), no_patches_per_type(a.no_patches_per_type), num_per_type(a.num_per_type), total_patches_per_type(a.total_patches_per_type), pot_starters(a.pot_starters)
{
    flat = a.flat;
    no_types = no_patches_per_type.getsize();
    
    int total_no_of_orients = 0;
    vector1<int> start_and_end_points(no_types + 1);
    start_and_end_points[0] = 0;

    for (int i = 0; i < no_patches_per_type.getsize(); i++)
    {
        total_no_of_orients += no_patches_per_type[i];
        start_and_end_points[i + 1] = total_no_of_orients;
    }



    int ht = no_types * (no_types+1) / 2;
    i1 = new int *[ht]; // allocate an array of 10 int pointers — these are our rows

    int hj = 0;
    int nn = 0;
    for (int i = 0; i < no_types; i++)
    {
        for (int j = i; j < no_types; j++)
        {


            pot_starters(i, j) = hj;

            i1[nn] = new int[1 + no_patches_per_type[i] * no_patches_per_type[j]];

            i1[nn][0] = no_patches_per_type[i] * no_patches_per_type[j];

            for (int k = 0; k < no_patches_per_type[i] * no_patches_per_type[j]; k++)
            {
                i1[nn][k + 1] = hj;
                hj++;
            }
            nn++;
        }
    }

    p = &i1[0];

    safe = false;

    int iter = 0;

    for (int i = 0; i < no_types; i++)
    {
        for (int j = i; j < no_types; j++)
        {
            // for(int i1  = 0 ; i1 < no_patches_per_type[i] ; i1++ ) {
            //     for(int j1 = 0 ; j1 < no_pathes_per_type[j] ; j1++) {

            int starti = start_and_end_points[i];
            int endi = start_and_end_points[i + 1];
            int startj = start_and_end_points[j];
            int endj = start_and_end_points[j + 1];

            for (int i2 = starti; i2 < endi; i2++)
            {
                for (int j1 = startj; j1 < endj; j1++)
                {
                    //mypot *pot1 = new mypot(orient(i2, 0), orient(i2, 1), orient(i2, 2), orient(j1, 0), orient(j1, 1), orient(j1, 2), params(iter, 0), params(iter, 1), params(iter, 2), 0.75);
                    potentialtheta3D *pot1;
                    if (flat)
                        pot1 = new KernFrenkelOnePatchFlatBottom(orient(i2, 0), orient(i2, 1), orient(i2, 2), orient(j1, 0), orient(j1, 1), orient(j1, 2), params(iter, 0), params(iter, 1), params(iter, 2), 0.5);
                    else
                        pot1 = new mypot(orient(i2, 0), orient(i2, 1), orient(i2, 2), orient(j1, 0), orient(j1, 1), orient(j1, 2), params(iter, 0), params(iter, 1), params(iter, 2), 0.2);

                    potential_bundle[iter] = pot1->clone();
                    delete pot1;
                    iter++;
                }
            }
        }
    }

}

inline int GeneralPatch::return_type(int i) {

    
    for(int j = 0 ; j < no_types; j++ ) 
        if(i < num_per_type[j] ) return j;


}

inline int GeneralPatch::return_type_patch(int i)
{

    // if (i < total_patches_per_type[0])
    //     return 0;

    for (int j = 0; j < no_types; j++)
        if (i < total_patches_per_type[j])
            return j ;
}

int GeneralPatch::num_patches(const int &i)
{
    int k1 = return_type(i);

    return no_patches_per_type[k1];
}

// inline int GeneralPatch::mapping_funcion_particles(int i, int j) {

// }

void GeneralPatch::UpdateIterator(const int &i, const int &j) {

int k1 = return_type(i);
int k2 = return_type(j);

// no_types*k1 - (k1*(1 + k1))/2. + k2

p = &(i1[ no_types*k1 - (k1*(1 + k1))/2 + k2 ]);
}


void GeneralPatch::UpdateIteratorSafe(const int &i, const int &j, int **q) {
    int k1 = return_type(i);
    int k2 = return_type(j);

    *q = (i1[no_types * k1 - (k1 * (1 + k1)) / 2 + k2]);
}

int GeneralPatch::get_total_patches(const int &N) {
    int num = 0 ;
    for(int i = 0  ; i < no_types; i++) {
    int k1 = num_per_type[i];
    int k2 = i > 0 ? num_per_type[i-1] : 0;
    
    num += no_patches_per_type[i] * (k1-k2);
    
    }
    return num;
 }

 void GeneralPatch::which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
 {
     int k1 = return_type(i);
     int k2 = return_type(j);

    int potstart = potn - pot_starters(k1,k2);

    int ori = no_patches_per_type[k2];

    int k3 = potstart / ori;
    int k4 = potstart % ori;

    int starti = k1 == 0 ? 0 : total_patches_per_type[k1 - 1];
    int startj = k2 == 0 ? 0 : total_patches_per_type[k2 - 1];

    int startpi = k1 == 0 ? 0 : num_per_type[k1 - 1];
    int startpj = k2 == 0 ? 0 : num_per_type[k2 - 1];



    wpi = starti + (i-startpi) * no_patches_per_type[k1] + k3;
    wpj = startj + (j-startpj) * no_patches_per_type[k2] + k4;

    //  wpi = i;
    //  wpj = j;
 }

 void GeneralPatch::which_particle(const int &wpi, const int &wpj, int &i, int &j)
 {

     int ki = return_type_patch(wpi);
     int kj = return_type_patch(wpj);

    
    int it = ki == 0 ? 0 : total_patches_per_type[ki - 1];
    int jt = kj == 0 ? 0 : total_patches_per_type[kj - 1];

    int starti = ki == 0 ? 0 : num_per_type[ki - 1];
    int startj = kj == 0 ? 0 : num_per_type[kj - 1];

    i = starti + (wpi - it) / no_patches_per_type[ki];
    j = startj + (wpj - jt) / no_patches_per_type[kj];
    //  i = wpi;
     //j = wpj;
 }

 int GeneralPatch::which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
 {
     int k1 = return_type(i);
     int k2 = return_type(j);

     int potstart = pot_starters(k1, k2);

     int starti = k1 == 0 ? 0 : total_patches_per_type[k1 - 1];
     int startj = k2 == 0 ? 0 : total_patches_per_type[k2 - 1];

     int patchi = (wpi - starti) % no_patches_per_type[k1];
     int patchj = (wpj - startj) % no_patches_per_type[k2];

     return potstart + patchi*no_patches_per_type[k2]+patchj;
 }

 void GeneralPatch::get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12) {
     d12 = potential_bundle[potn]->interaction_distance;
     ang12 = params(potn, 2);

     int k1 = return_type(i);
     int k2 = return_type(j);

     int potstart = potn - pot_starters(k1, k2);

     int ori = no_patches_per_type[k2];



     int k3 = potstart / ori;
     int k4 = potstart % ori;


     int starti = 0;
     int startj = 0;
     for(int i1 = 0 ; i1 < k1 ; i1++)
        starti += no_patches_per_type[i1];

    for (int j1 = 0; j1 < k2; j1++)
         startj += no_patches_per_type[j1];


    nxb1 = orient(starti + k3, 0);
    nyb1 = orient(starti + k3, 1);
    nzb1 = orient(starti + k3, 2);

    nxb2 = orient(startj + k4, 0);
    nyb2 = orient(startj + k4, 1);
    nzb2 = orient(startj + k4, 2);
 }

#endif /* GENERALPATCH_CPP */
