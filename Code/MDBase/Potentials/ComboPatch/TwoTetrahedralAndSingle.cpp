#ifndef TWOTETRAHEDRALANDSINGLE_CPP
#define TWOTETRAHEDRALANDSINGLE_CPP

TwoTetrahedralAndSingle::TwoTetrahedralAndSingle(matrix<double> &paramss, int ntt, int nss, int nff) : ComboPatch(57), params(paramss), v(matrix<double>(4, 3))
{

    if (params.getnrows() != 6)
        error("number of rows in parameter matrix is incorrect");
    if (params.getncols() != 3)
        error("number of cols in parameter matrix is incorrect");
    nt = ntt;
    ns = nss;
    nf = nff;
    i1 = new int[17];
    i2 = new int[17];
    i3 = new int[5];
    i4 = new int[17];
    i5 = new int[5];
    i6 = new int[2];

    i1[0] = 16;
    i2[0] = 16;
    i3[0] = 4;
    i4[0] = 16;
    i5[0] = 4;
    i6[0] = 1;

    int iter = 0;
    for (int i = 0; i < 16; i++)
        i1[i + 1] = iter + i;

    iter += 16;

    for (int i = 0; i < 16; i++)
        i2[i + 1] = iter + i;

    iter += 16;

    for (int i = 0; i < 4; i++)
        i3[i + 1] = iter + i;

    iter += 4;

    for (int i = 0; i < 16; i++)
        i4[i + 1] = iter + i;

    iter += 16;

    for (int i = 0; i < 4; i++)
        i5[i + 1] = iter + i;

    iter += 4;

    for (int i = 0; i < 1; i++)
        i6[i + 1] = iter + i;

    p = &i1;

    safe = false;

    //matrix<double> v(4, 3);
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

    iter = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), params(0, 0), params(0, 1), params(0, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), params(1, 0), params(1, 1), params(1, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), params(2, 0), params(2, 1), params(2, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), params(3, 0), params(3, 1), params(3, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), params(4, 0), params(4, 1), params(4, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), params(5, 0), params(5, 1), params(5, 2), 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
}


int TwoTetrahedralAndSingle::num_patches(const int &i)
{
    if (i < ns)
        return 4;
    else
        return 1;
}

inline int TwoTetrahedralAndSingle::mapping_funcion_particles(int i, int j)
{
    if (i < nt)
    {
        if (j < nt)
        {
            return 1;
        }
        else if (j < ns)
        {
            return 2;
        }
        else
        {
            return 3;
        }
    }
    else if (i < ns)
    {
        if (j < nt)
        {
            return 2;
        }
        else if (j < ns)
        {
            return 4;
        }
        else
        {
            return 5;
        }
    }
    else
    {
        if (j < nt)
        {
            return 3;
        }
        else if (j < ns)
        {
            return 5;
        }
        else
        {
            return 6;
        }
    }
}
void TwoTetrahedralAndSingle::UpdateIterator(const int &i, const int &j)
{
    int we_are = mapping_funcion_particles(i, j);

    switch (we_are)
    {
    case 1:
        p = &i1;
        break;

    case 2:
        p = &i2;
        break;

    case 3:
        p = &i3;
        break;

    case 4:
        p = &i4;
        break;

    case 5:
        p = &i5;
        break;

    case 6:
        p = &i6;
        break;
    }
}
void TwoTetrahedralAndSingle::UpdateIteratorSafe(const int &i, const int &j, int **q)
{
    int we_are = mapping_funcion_particles(i, j);

    switch (we_are)
    {
    case 1:
        *q = i1;
        break;

    case 2:
        *q = i2;
        break;

    case 3:
        *q = i3;
        break;

    case 4:
        *q = i4;
        break;

    case 5:
        *q = i5;
        break;

    case 6:
        *q = i6;
        break;
    }
}

int TwoTetrahedralAndSingle::get_total_patches(const int &N) { return 4 * nt + 4 * (ns - nt) + nf - ns - nt; }

void TwoTetrahedralAndSingle::which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
{
    int we_are = mapping_funcion_particles(i, j);

    switch (we_are)
    {
    case 1:
    {
        int k1 = potn / 4;
        int k2 = potn % 4;
        wpi = i * 4 + k1;
        wpj = j * 4 + k2;
        break;
    }
    case 2:
    {
        int potn2 = potn - 16;
        int k1 = potn2 / 4;
        int k2 = potn2 % 4;
        wpi = i * 4 + k1;
        wpj = j * 4 + k2;
        break;
    }
    case 3:
    {
        int k1 = (potn - 32) % 4;
        if (i < j)
        {
            wpi = i * 4 + k1;
            wpj = ns * 4 + (j - ns);
        }
        else
        {
            wpj = j * 4 + k1;
            wpi = ns * 4 + (i - ns);
        }
        break;
    }
    case 4:
    {
        int potn2 = potn - 36;
        int k1 = potn2 / 4;
        int k2 = potn2 % 4;
        wpi = i * 4 + k1;
        wpj = j * 4 + k2;
        break;
    }
    case 5:
    {
        int k1 = (potn - 52) % 4;
        if (i < j)
        {
            wpi = i * 4 + k1;
            wpj = ns * 4 + (j - ns);
        }
        else
        {
            wpj = j * 4 + k1;
            wpi = ns * 4 + (i - ns);
        }
        break;
    }
    case 6:
    {
        wpi = ns * 4 + (i - ns);
        wpj = ns * 4 + (j - ns);
        break;
    }
    }
}

void TwoTetrahedralAndSingle::which_particle(const int &wpi, const int &wpj, int &i, int &j)
{
    if (wpi < 4 * ns && wpj < 4 * ns)
    {
        i = wpi / 4;
        j = wpj / 4;
    }
    else if (wpi >= 4 * ns && wpj >= 4 * ns)
    {
        i = wpi - 4 * ns + ns;
        j = wpj - 4 * ns + ns;
    }
    else if (wpi < 4 * ns && wpj >= 4 * ns)
    {
        i = wpi / 4;
        j = wpj - 4 * ns + ns;
    }
    else
    {
        i = wpi - 4 * ns + ns;
        j = wpj / 4;
    }
}
int TwoTetrahedralAndSingle::which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
{
    int we_are = mapping_funcion_particles(i, j);
    switch (we_are)
    {
    case 1:
        return (wpi % 4) * 4 + (wpj % 4);
        break;

    case 2:
        return 16 + (wpi % 4) * 4 + (wpj % 4);
        break;

    case 3:
        if (i < j)
            return 32 + wpi % 4;
        else
            return 32 + wpj % 4;
        break;

    case 4:
        return 36 + (wpi % 4) * 4 + (wpj % 4);
        break;

    case 5:
        if (i < j)
            return 52 + wpi % 4;
        else
            return 52 + wpj % 4;
        break;

    case 6:
        return 56;
        break;
        
    default:
        error("mapping function failed in which potential");
        return 0;
        break;
    }

  
}

void TwoTetrahedralAndSingle::get_params(const int &i, const int &j, const int &potn2, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12)
{

    int we_are = mapping_funcion_particles(i, j);
    int i1;
    int i2;
    int potn;

    switch (we_are)
    {
    case 1:
        i1 = potn2 / 4;
        i2 = potn2 % 4;
        d12 = potential_bundle[potn2]->interaction_distance;
        ang12 = params(0, 2);
        nxb1 = v(i1, 0);
        nyb1 = v(i1, 1);
        nzb1 = v(i1, 2);
        nxb2 = v(i2, 0);
        nyb2 = v(i2, 1);
        nzb2 = v(i2, 2);
        break;

    case 2:
        potn = potn2 - 16;
        i1 = potn / 4;
        i2 = potn % 4;
        d12 = potential_bundle[potn2]->interaction_distance;
        ang12 = params(1, 2);
        nxb1 = v(i1, 0);
        nyb1 = v(i1, 1);
        nzb1 = v(i1, 2);
        nxb2 = v(i2, 0);
        nyb2 = v(i2, 1);
        nzb2 = v(i2, 2);
        break;

    case 3:
        potn = potn2 - 32;
        if (i < j)
        {
            i1 = potn % 4;
            i2 = 0;
        }
        else
        {
            i2 = potn % 4;
            i1 = 0;
        }
        d12 = potential_bundle[potn2]->interaction_distance;
        ;
        ang12 = params(2, 2);
        nxb1 = v(i1, 0);
        nyb1 = v(i1, 1);
        nzb1 = v(i1, 2);
        nxb2 = v(i2, 0);
        nyb2 = v(i2, 1);
        nzb2 = v(i2, 2);
        break;

    case 4:
        potn = potn2 - 36;
        i1 = potn / 4;
        i2 = potn % 4;
        d12 = potential_bundle[potn2]->interaction_distance;
        ang12 = params(3, 2);
        nxb1 = v(i1, 0);
        nyb1 = v(i1, 1);
        nzb1 = v(i1, 2);
        nxb2 = v(i2, 0);
        nyb2 = v(i2, 1);
        nzb2 = v(i2, 2);
        break;

    case 5:
        potn = potn2 - 52;
        if (i < j)
        {
            i1 = potn % 4;
            i2 = 0;
        }
        else
        {
            i2 = potn % 4;
            i1 = 0;
        }
        d12 = potential_bundle[potn2]->interaction_distance;
        ang12 = params(4, 2);
        nxb1 = v(i1, 0);
        nyb1 = v(i1, 1);
        nzb1 = v(i1, 2);
        nxb2 = v(i2, 0);
        nyb2 = v(i2, 1);
        nzb2 = v(i2, 2);
        break;

    case 6:
        potn = potn2 - 56;
        i1 = 0;
        i2 = 0;
        d12 = potential_bundle[potn2]->interaction_distance;
        ang12 = params(5, 2);
        nxb1 = v(i1, 0);
        nyb1 = v(i1, 1);
        nzb1 = v(i1, 2);
        nxb2 = v(i2, 0);
        nyb2 = v(i2, 1);
        nzb2 = v(i2, 2);
        break;
    }
}

#endif /* TWOTETRAHEDRALANDSINGLE_CPP */
