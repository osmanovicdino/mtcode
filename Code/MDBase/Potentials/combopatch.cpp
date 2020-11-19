#ifndef COMBOPATCH_CPP
#define COMBOPATCH_CPP

SingPatch::SingPatch(double strr, double disss, double angg) : ComboPatch(1), ang(angg), dis(disss), str(strr)
{
    i1 = new int[2];
    i1[0] = 1;
    i1[1] = 0;
    p = &i1;


    KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(nx, ny, nz, nx, ny, nz, str, dis, ang, 0.75);

    potential_bundle[0] = pot1->clone();


    delete pot1;
    //potential_bundle[0] = pot1;
    //vector1<potentialtheta3D *> pots(1);
}

TetrahedralPatch::TetrahedralPatch(double strr, double disss, double angg) : ComboPatch(16), ang(angg), dis(disss), str(strr) {
    i1 = new int[17];
    i1[0] = 16;
    for(int i = 0 ; i < 16 ; i++)
        i1[i+1] = i;
    p = &i1;

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
            KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), str, dis, ang, 0.75);

            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }


}

void TetrahedralPatch::get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12)
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
        d12 = dis;
        ang12 = ang;
    
}

TetrahedralWithSingle::TetrahedralWithSingle(double strrtt, double disstt, double anggtt, double strrts, double dissts, double anggts, double strrss, double dissss, double anggss, int ntt, int nss) : ComboPatch(21) 
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
    for(int i = 0  ; i < 4 ; i++)
        i2[i+1] = 16+i;
   
        i3[1] = 20;
    p = &i1;

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
            KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), strtt, distt, angtt, 0.75);

            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), strts, dists, angts, 0.75);

            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            KernFrenkelOnePatch2 *pot1 = new KernFrenkelOnePatch2(v(i, 0), v(i, 1), v(i, 2), v(j, 0), v(j, 1), v(j, 2), strss, disss, angss, 0.75);

            potential_bundle[iter] = pot1->clone();
            delete pot1;
            iter++;
        }
    }

}

void TetrahedralWithSingle::get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12) {
    if(i < nt && j < nt ) {
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
        d12 = distt;
        ang12 = angtt;
    }
    else if(i >= nt && j >= nt) {
        nxb1 = nx1;
        nyb1 = ny1;
        nzb1 = nz1;
        nxb2 = nx1;
        nyb2 = ny1;
        nzb2 = nz1;
        d12 =  disss;
        ang12 = angss;
    }
    else if(i<nt && j >= nt) {
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
        d12 = dists;
        ang12 = angts;

    }
    else {
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

#endif /* COMBOPATCH_CPP */
