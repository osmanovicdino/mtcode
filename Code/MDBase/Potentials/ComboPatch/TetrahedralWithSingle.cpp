#ifndef TETRAHEDRALWITHSINGLE_CPP
#define TETRAHEDRALWITHSINGLE_CPP

TetrahedralWithSingle::TetrahedralWithSingle(double strrtt, double disstt, double anggtt, double strrts, double dissts, double anggts, double strrss, double dissss, double anggss, int ntt, int nss) : ComboPatch(21), v(matrix<double>(4, 3))
{
    nt = ntt;
    ns = nss;
    i1 = new int[17];
    i2 = new int[5];
    i3 = new int[2];
    i1[0] = 16;
    i2[0] = 4;
    i3[0] = 1;
    for (int i = 0; i < 16; i++)
        i1[i + 1] = i;
    for (int i = 0; i < 4; i++)
        i2[i + 1] = 16 + i;

    i3[1] = 20;
    p = &i1;

    safe = false;

    strtt = strrtt;
    strts = strrts;
    strss = strrss;

    distt = disstt;
    dists = dissts;
    disss = dissss;

    angtt = anggtt;
    angts = anggts;
    angss = anggss;

    matrix<double> v(4, 3);
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

    int iter = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), strtt, distt, angtt, 0.75);

            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), strts, dists, angts, 0.75);

            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), strss, disss, angss, 0.75);

            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
}


int TetrahedralWithSingle::num_patches(const int &i)
{
    if (i < nt)
    {
        return 4;
    }
    else if (i >= nt)
    {
        return 1;
    }
    else
    {
        error("out of bounds");
    }
}

void TetrahedralWithSingle::UpdateIterator(const int &i, const int &j)
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

void TetrahedralWithSingle::UpdateIteratorSafe(const int &i, const int &j, int **q)
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

int TetrahedralWithSingle::get_total_patches(const int &N) { return 4 * nt + ns; }

void TetrahedralWithSingle::which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
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
        wpi = nt * 4 + (i - nt);
        wpj = nt * 4 + (j - nt);
    }
    else if (i < nt && j >= nt)
    {
        int k1 = potn % 4;
        wpi = i * 4 + k1;
        wpj = nt * 4 + (j - nt);
    }
    else
    {
        int k2 = potn % 4;
        wpi = nt * 4 + (i - nt);
        wpj = j * 4 + k2;
    }
}

void TetrahedralWithSingle::which_particle(const int &wpi, const int &wpj, int &i, int &j)
{
    if (wpi < 4 * nt && wpj < 4 * nt)
    {
        i = wpi / 4;
        j = wpj / 4;
    }
    else if (wpi >= 4 * nt && wpj >= 4 * nt)
    {
        i = wpi - 4 * nt + nt;
        j = wpj - 4 * nt + nt;
    }
    else if (wpi < 4 * nt && wpj >= 4 * nt)
    {
        i = wpi / 4;
        j = wpj - 4 * nt + nt;
    }
    else
    {
        i = wpi - 4 * nt + nt;
        j = wpj / 4;
    }
}

int TetrahedralWithSingle::which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
{
    if (wpi < 4 * nt && wpj < 4 * nt)
    {
        return (wpi % 4) * 4 + (wpj % 4);
    }
    else if (wpi >= 4 * nt && wpj >= 4 * nt)
    {
        return 20;
    }
    else if (wpi < 4 * nt && wpj >= 4 * nt)
    {
        return 16 + wpi % 4;
    }
    else
    {
        return 16 + wpj % 4;
    }
}

void TetrahedralWithSingle::get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12)
{
    if (i < nt && j < nt)
    {
        if (potn == 0)
        {
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
        }
        else if (potn == 1)
        {
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
            nxb2 = nx2;
            nyb2 = ny2;
            nzb2 = nz2;
        }
        else if (potn == 2)
        {
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
            nxb2 = nx3;
            nyb2 = ny3;
            nzb2 = nz3;
        }
        else if (potn == 3)
        {
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
            nxb2 = nx4;
            nyb2 = ny4;
            nzb2 = nz4;
        }
        else if (potn == 4)
        {
            nxb1 = nx2;
            nyb1 = ny2;
            nzb1 = nz2;
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
        }
        else if (potn == 5)
        {
            nxb1 = nx2;
            nyb1 = ny2;
            nzb1 = nz2;
            nxb2 = nx2;
            nyb2 = ny2;
            nzb2 = nz2;
        }
        else if (potn == 6)
        {
            nxb1 = nx2;
            nyb1 = ny2;
            nzb1 = nz2;
            nxb2 = nx3;
            nyb2 = ny3;
            nzb2 = nz3;
        }
        else if (potn == 7)
        {
            nxb1 = nx2;
            nyb1 = ny2;
            nzb1 = nz2;
            nxb2 = nx4;
            nyb2 = ny4;
            nzb2 = nz4;
        }
        else if (potn == 8)
        {
            nxb1 = nx3;
            nyb1 = ny3;
            nzb1 = nz3;
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
        }
        else if (potn == 9)
        {
            nxb1 = nx3;
            nyb1 = ny3;
            nzb1 = nz3;
            nxb2 = nx2;
            nyb2 = ny2;
            nzb2 = nz2;
        }
        else if (potn == 10)
        {
            nxb1 = nx3;
            nyb1 = ny3;
            nzb1 = nz3;
            nxb2 = nx3;
            nyb2 = ny3;
            nzb2 = nz3;
        }
        else if (potn == 11)
        {
            nxb1 = nx3;
            nyb1 = ny3;
            nzb1 = nz3;
            nxb2 = nx4;
            nyb2 = ny4;
            nzb2 = nz4;
        }
        else if (potn == 12)
        {
            nxb1 = nx4;
            nyb1 = ny4;
            nzb1 = nz4;
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
        }
        else if (potn == 13)
        {
            nxb1 = nx4;
            nyb1 = ny4;
            nzb1 = nz4;
            nxb2 = nx4;
            nyb2 = ny2;
            nzb2 = nz2;
        }
        else if (potn == 14)
        {
            nxb1 = nx4;
            nyb1 = ny4;
            nzb1 = nz4;
            nxb2 = nx3;
            nyb2 = ny3;
            nzb2 = nz3;
        }
        else if (potn == 15)
        {
            nxb1 = nx4;
            nyb1 = ny4;
            nzb1 = nz4;
            nxb2 = nx4;
            nyb2 = ny4;
            nzb2 = nz4;
        }
        else
        {
            cout << i << " " << j << endl;
            cout << "potential number: " << potn << endl;
            error("potn out of bounds in tetrahedral get_params");
        }
        d12 = potential_bundle[potn]->interaction_distance;
        ang12 = angtt;
    }
    else if (i >= nt && j >= nt)
    {
        nxb1 = nx1;
        nyb1 = ny1;
        nzb1 = nz1;
        nxb2 = nx1;
        nyb2 = ny1;
        nzb2 = nz1;
        d12 = disss;
        ang12 = angss;
    }
    else if (i < nt && j >= nt)
    {
        if (potn == 16)
        {
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
        }
        else if (potn == 17)
        {
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
            nxb1 = nx2;
            nyb1 = ny2;
            nzb1 = nz2;
        }
        else if (potn == 18)
        {
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
            nxb1 = nx3;
            nyb1 = ny3;
            nzb1 = nz3;
        }
        else if (potn == 19)
        {
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
            nxb1 = nx4;
            nyb1 = ny4;
            nzb1 = nz4;
        }
        d12 = potential_bundle[potn]->interaction_distance;
        ang12 = angts;
    }
    else
    {
        if (potn == 16)
        {
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
            nxb2 = nx1;
            nyb2 = ny1;
            nzb2 = nz1;
        }
        else if (potn == 17)
        {
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
            nxb2 = nx2;
            nyb2 = ny2;
            nzb2 = nz2;
        }
        else if (potn == 18)
        {
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
            nxb2 = nx3;
            nyb2 = ny3;
            nzb2 = nz3;
        }
        else if (potn == 19)
        {
            nxb1 = nx1;
            nyb1 = ny1;
            nzb1 = nz1;
            nxb2 = nx4;
            nyb2 = ny4;
            nzb2 = nz4;
        }
        d12 = dists;
        ang12 = angts;
    }
}

#endif /* TETRAHEDRALWITHSINGLE_CPP */
