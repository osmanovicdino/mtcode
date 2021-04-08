#ifndef BIVALENTPATCH_CPP
#define BIVALENTPATCH_CPP

BivalentPatch::BivalentPatch(double strr, double disss, double angg) : ComboPatch(4), ang(angg), dis(disss), str(strr)
{
    i1 = new int[5];
    i1[0] = 4;
    for (int i = 0; i < 4; i++)
        i1[i + 1] = i;
    p = &i1;

    
    matrix<double> v(2,3);
    v(0, 0) = nx1;
    v(0, 1) = ny1;
    v(0, 2) = nz1;

    //vector1<double> v2(3);

    v(1, 0) = nx2;
    v(1, 1) = ny2;
    v(1, 2) = nz2;

    int iter = 0;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), str, dis, ang, 0.75);
            if(pot1->interaction_distance > max_check) max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
}

BivalentPatch::BivalentPatch(const BivalentPatch &old) : ComboPatch(4), ang(old.ang), dis(old.dis), str(old.str)
{
    i1 = new int[5];
    i1[0] = 4;
    for (int i = 0; i < 4; i++)
        i1[i + 1] = i;
    p = &i1;

    matrix<double> v(2, 3);
    v(0, 0) = nx1;
    v(0, 1) = ny1;
    v(0, 2) = nz1;

    //vector1<double> v2(3);

    v(1, 0) = nx2;
    v(1, 1) = ny2;
    v(1, 2) = nz2;

    int iter = 0;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), str, dis, ang, 0.75);

            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
}

int BivalentPatch::num_patches(const int &i)
{
    return 2;
}
void BivalentPatch::UpdateIterator(const int &i, const int &j)
{
    //do nothing (only one patch)
}
void BivalentPatch::UpdateIteratorSafe(const int &i, const int &j, int **q)
{
    *q = i1;
}

int BivalentPatch::get_total_patches(const int &N) { return 2 * N; }
void BivalentPatch::which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
{
    int k1 = potn / 2;
    int k2 = potn % 2;
    // i*4 + k1;
    // j*4 + k2;
    wpi = i * 2 + k1;
    wpj = j * 2 + k2;
}

void BivalentPatch::which_particle(const int &wpi, const int &wpj, int &i, int &j)
{
    i = wpi / 2;
    j = wpj / 2;
}
int BivalentPatch::which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
{
    return (wpi % 2) * 2 + (wpj % 2);
}

inline void BivalentPatch::get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12)
{
    d12 = potential_bundle[potn]->interaction_distance;
    ang12 = ang;

    switch (potn)
    {
    case 0:
        nxb1 = nx1;
        nyb1 = ny1;
        nzb1 = nz1;
        nxb2 = nx1;
        nyb2 = ny1;
        nzb2 = nz1;
        break;

    case 1:
        nxb1 = nx1;
        nyb1 = ny1;
        nzb1 = nz1;
        nxb2 = nx2;
        nyb2 = ny2;
        nzb2 = nz2;
        break;

    case 2:
        nxb1 = nx2;
        nyb1 = ny2;
        nzb1 = nz2;
        nxb2 = nx1;
        nyb2 = ny1;
        nzb2 = nz1;
        break;

    case 3:
        nxb1 = nx2;
        nyb1 = ny2;
        nzb1 = nz2;
        nxb2 = nx2;
        nyb2 = ny2;
        nzb2 = nz2;
        break;
    default:
        cout << i << " " << j << endl;
        cout << "potential number: " << potn << endl;
        error("potn out of bounds in tetrahedral get_params");
        break;
    }
}

#endif /* BIVALENTPATCH_CPP */
