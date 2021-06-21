#ifndef TETRAHEDRALWITHBIVALENT_CPP
#define TETRAHEDRALWITHBIVALENT_CPP

TetrahedralWithBivalent::TetrahedralWithBivalent(matrix<double> &params, int ntt, int nbb) : ComboPatch(28), v(matrix<double>(4, 3)), v2(matrix<double>(2,3)),params2(params)
{
    if(params.getNsafe()!=28) error("param list for Tetrahedron with Bivalent not of correct dimension (should be 28)");
    nt = ntt;
    nb = nbb;
    i1 = new int[17];
    i2 = new int[9];
    i3 = new int[5];
    i1[0] = 16;
    i2[0] = 8;
    i3[0] = 4;
    for (int i = 0; i < 16; i++)
        i1[i + 1] = i;
    for (int i = 0; i < 8; i++)
        i2[i + 1] = 16 + i;
    for (int i = 0; i < 4; i++)
        i3[i + 1] = 24 + i;

    p = &i1;

    safe = false;

    v(0, 0) = nx1;
    v(0, 1) = ny1;
    v(0, 2) = nz1;

    //vector1<double> v2(3);

    v(1, 0) = nx2;
    v(1, 1) = ny2;
    v(1, 2) = nz2;

    //vector1<double> v3(3);

    v(2, 0) = nx3;
    v(2, 1) = ny3;
    v(2, 2) = nz3;

    //vector1<double> v4(3);

    v(3, 0) = nx4;
    v(3, 1) = ny4;
    v(3, 2) = nz4;

    v2(0, 0) = nx5;
    v2(0, 1) = ny5;
    v2(0, 2) = nz5;

    v2(1, 0) = nx6;
    v2(1, 1) = ny6;
    v2(1, 2) = nz6;

    int iter = 0;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), params(iter, 0), params(iter, 1), params(iter, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v2(j, 0), v2(j, 1), v2(j, 2), params(iter, 0), params(iter, 1), params(iter, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            mypot *pot1 = new mypot(v2(i, 0), v2(i, 1), v2(i, 2), v2(j, 0), v2(j, 1), v2(j, 2), params(iter, 0), params(iter, 1), params(iter, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
}

int TetrahedralWithBivalent::num_patches(const int &i)
{
    if (i < nt)
    {
        return 4;
    }
    else if (i >= nt)
    {
        return 2;
    }
    else
    {
        error("out of bounds");
        return 0;
    }
}

void TetrahedralWithBivalent::UpdateIterator(const int &i, const int &j)
{
    if (i < nt && j < nt)
    {
        p = &i1;
    }
    else if (i >= nt && j >= nt)
    {
        p = &i3;
    }
    else
    {
        p = &i2;
    }
}

void TetrahedralWithBivalent::UpdateIteratorSafe(const int &i, const int &j, int **q)
{
    if (i < nt && j < nt)
    {
        *q = i1;
    }
    else if (i >= nt && j >= nt)
    {
        *q = i3;
    }
    else
    {
        *q = i2;
    }
}

int TetrahedralWithBivalent::get_total_patches(const int &N) { return 4 * nt + 2*(nb-nt); }

void TetrahedralWithBivalent::which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
{
    if (i < nt && j < nt)
    {
        int k1 = potn / 4;
        int k2 = potn % 4;
        // i*4 + k1;
        // j*4 + k2;
        wpi = i * 4 + k1;
        wpj = j * 4 + k2;
    }
    else if (i >= nt && j >= nt)
    {
        int k1 = (potn-24)/2;
        int k2 = (potn-24) %2;

        wpi = nt * 4 + (i - nt)*2 + k1;
        wpj = nt * 4 + (j - nt)*2 + k2;
    }
    else if (i < nt && j >= nt)
    {
        int k1 = (potn-16) / 2;
        int k2 = (potn-16) % 2;
        wpi = i * 4 + k1;
        wpj = nt * 4 + (j - nt)*2 + k2;
    }
    else
    {
        int k1 = (potn - 16) / 2;
        int k2 = (potn - 16) % 2;
        wpi = nt * 4 + (i - nt)*2 + k2;
        wpj = j * 4 + k1;
    }
}

void TetrahedralWithBivalent::which_particle(const int &wpi, const int &wpj, int &i, int &j)
{
    if (wpi < 4 * nt && wpj < 4 * nt)
    {
        i = wpi / 4;
        j = wpj / 4;
    }
    else if (wpi >= 4 * nt && wpj >= 4 * nt)
    {
        i = (wpi - 4 * nt)/2 + nt;
        j = (wpj - 4 * nt)/2 + nt;
    }
    else if (wpi < 4 * nt && wpj >= 4 * nt)
    {
        i = wpi / 4;
        j = (wpj - 4 * nt)/2 + nt;
    }
    else
    {
        i = (wpi - 4 * nt)/2 + nt;
        j = wpj / 4;
    }
}

int TetrahedralWithBivalent::which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
{
    if (wpi < 4 * nt && wpj < 4 * nt)
    {
        return (wpi % 4) * 4 + (wpj % 4);
    }
    else if (wpi >= 4 * nt && wpj >= 4 * nt)
    {
        return 24 + ((wpi-4*nt) % 2) *2 + ((wpj-4*nt) % 2);
    }
    else if (wpi < 4 * nt && wpj >= 4 * nt)
    {
        return 16 + (wpi % 4) * 2 + ((wpj - 4 * nt) % 2);
    }
    else
    {
        return 16 + (wpj % 4) * 2 + ((wpi - 4 * nt) % 2);
    }
}

void TetrahedralWithBivalent::get_params(const int &i1, const int &j1, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12)
{
    d12 = potential_bundle[potn]->interaction_distance;
    ang12 =params2(potn,2);

    if(potn < 16) {
        int i = potn/4;
        int j = potn%4;
        nxb1 = v(i,0);
        nyb1 = v(i,1);
        nzb1 = v(i,2);
        nxb2 = v(j,0);
        nyb2 = v(j,1);
        nzb2 = v(j,2);
    }
    else if(potn < 24) {
        int i = (potn - 16 ) /2;
        int j = (potn - 16 ) %2;
        nxb1 = v(i, 0);
        nyb1 = v(i, 1);
        nzb1 = v(i, 2);
        nxb2 = v2(j, 0);
        nyb2 = v2(j, 1);
        nzb2 = v2(j, 2);
    }
    else{
        int i = (potn - 24) / 2;
        int j = (potn - 24) % 2;
        nxb1 = v2(i, 0);
        nyb1 = v2(i, 1);
        nzb1 = v2(i, 2);
        nxb2 = v2(j, 0);
        nyb2 = v2(j, 1);
        nzb2 = v2(j, 2);
    }
}

#endif /* TETRAHEDRALWITHBIVALENT_CPP */