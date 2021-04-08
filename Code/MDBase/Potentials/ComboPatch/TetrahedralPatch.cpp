#ifndef TETRAHEDRALPATCH_CPP
#define TETRAHEDRALPATCH_CPP

TetrahedralPatch::TetrahedralPatch(double strr, double disss, double angg) : ComboPatch(16), ang(angg), dis(disss), str(strr), v(matrix<double>(4, 3))
{
    i1 = new int[17];
    i1[0] = 16;
    for (int i = 0; i < 16; i++)
        i1[i + 1] = i;
    p = &i1;

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

    int iter = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mypot *pot1 = new mypot(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), str, dis, ang, 0.75);
            if (pot1->interaction_distance > max_check)
                max_check = pot1->interaction_distance;
            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
}


int TetrahedralPatch::num_patches(const int &i)
{
    return 4;
}
void TetrahedralPatch::UpdateIterator(const int &i, const int &j)
{
    //do nothing (only one patch)
}
void TetrahedralPatch::UpdateIteratorSafe(const int &i, const int &j, int **q)
{
    *q = i1;
}

int TetrahedralPatch::get_total_patches(const int &N) { return 4 * N; }
void TetrahedralPatch::which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
{
    int k1 = potn / 4;
    int k2 = potn % 4;
    // i*4 + k1;
    // j*4 + k2;
    wpi = i * 4 + k1;
    wpj = j * 4 + k2;
}
void TetrahedralPatch::which_particle(const int &wpi, const int &wpj, int &i, int &j)
{
    i = wpi / 4;
    j = wpj / 4;
}
int TetrahedralPatch::which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
{
    return (wpi % 4) * 4 + (wpj % 4);
}
// void TetrahedralPatch::get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12);

inline void TetrahedralPatch::get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12)
{
    // int i1= potn/4;
    // int i2 =potn%4;
    // nxb1 = v(i1,0);
    // nyb1 = v(i1, 1);
    // nzb1 = v(i1, 2);
    // nxb2 = v(i2, 0);
    // nyb2 = v(i2, 1);
    // nzb2 = v(i2, 2);
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
        nxb1 = nx1;
        nyb1 = ny1;
        nzb1 = nz1;
        nxb2 = nx3;
        nyb2 = ny3;
        nzb2 = nz3;
        break;

    case 3:
        nxb1 = nx1;
        nyb1 = ny1;
        nzb1 = nz1;
        nxb2 = nx4;
        nyb2 = ny4;
        nzb2 = nz4;
        break;

    case 4:
        nxb1 = nx2;
        nyb1 = ny2;
        nzb1 = nz2;
        nxb2 = nx1;
        nyb2 = ny1;
        nzb2 = nz1;
        break;

    case 5:
        nxb1 = nx2;
        nyb1 = ny2;
        nzb1 = nz2;
        nxb2 = nx2;
        nyb2 = ny2;
        nzb2 = nz2;
        break;

    case 6:
        nxb1 = nx2;
        nyb1 = ny2;
        nzb1 = nz2;
        nxb2 = nx3;
        nyb2 = ny3;
        nzb2 = nz3;
        break;

    case 7:
        nxb1 = nx2;
        nyb1 = ny2;
        nzb1 = nz2;
        nxb2 = nx4;
        nyb2 = ny4;
        nzb2 = nz4;
        break;

    case 8:
        nxb1 = nx3;
        nyb1 = ny3;
        nzb1 = nz3;
        nxb2 = nx1;
        nyb2 = ny1;
        nzb2 = nz1;
        break;

    case 9:
        nxb1 = nx3;
        nyb1 = ny3;
        nzb1 = nz3;
        nxb2 = nx2;
        nyb2 = ny2;
        nzb2 = nz2;
        break;

    case 10:
        nxb1 = nx3;
        nyb1 = ny3;
        nzb1 = nz3;
        nxb2 = nx3;
        nyb2 = ny3;
        nzb2 = nz3;
        break;

    case 11:
        nxb1 = nx3;
        nyb1 = ny3;
        nzb1 = nz3;
        nxb2 = nx4;
        nyb2 = ny4;
        nzb2 = nz4;
        break;

    case 12:
        nxb1 = nx4;
        nyb1 = ny4;
        nzb1 = nz4;
        nxb2 = nx1;
        nyb2 = ny1;
        nzb2 = nz1;
        break;

    case 13:
        nxb1 = nx4;
        nyb1 = ny4;
        nzb1 = nz4;
        nxb2 = nx4;
        nyb2 = ny2;
        nzb2 = nz2;
        break;

    case 14:
        nxb1 = nx4;
        nyb1 = ny4;
        nzb1 = nz4;
        nxb2 = nx3;
        nyb2 = ny3;
        nzb2 = nz3;
        break;

    case 15:
        nxb1 = nx4;
        nyb1 = ny4;
        nzb1 = nz4;
        nxb2 = nx4;
        nyb2 = ny4;
        nzb2 = nz4;
        break;

    default:
        cout << i << " " << j << endl;
        cout << "potential number: " << potn << endl;
        error("potn out of bounds in tetrahedral get_params");
        break;
    }
}

#endif /* TETRAHEDRALPATCH_CPP */
