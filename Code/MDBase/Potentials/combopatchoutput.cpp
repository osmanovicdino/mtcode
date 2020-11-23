#ifndef COMBOPATCHOUTPUT_CPP
#define COMBOPATCHOUTPUT_CPP

void ComboPatch::CreateFiles() {
    matrix<double> v(1, 3);
    v(0, 0) = 1.0;
    v(0, 1) = 0.0;
    v(0, 2) = 0.0;

    outfunc(v, "ori");

    matrix<double> col(1, 3);
    col(0, 0) = 0.87;
    col(0, 1) = 0.94;
    col(0, 2) = 1.;

    outfunc(col, "col");
}

void TetrahedralPatch::CreateFiles() {
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

    outfunc(v,"ori");

    matrix<double> col(1, 3);
    col(0, 0) = 0.87;
    col(0, 1) = 0.94;
    col(0, 2) = 1.;

    outfunc(col, "col");
}

void TetrahedralWithSingle::CreateFiles()
{
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

    outfunc(v, "ori");

    matrix<double> v2(1,3);

    v2(0, 0) = 1.0;
    v2(0, 1) = 0.0;
    v2(0, 2) = 0.0;

    outfunc(v2, "ori2");

    matrix<double> col(nt+ns, 3);
    for(int i = 0 ; i < nt ; i++)
    {
    col(i, 0) = 0.87;
    col(i, 1) = 0.94;
    col(i, 2) = 1.;
    }
    for(int i = nt  ; i < nt+ns ; i++) {
    col(i, 0) = 1.;
    col(i, 1) = 0.5;
    col(i, 2) = 0.;
    }

    outfunc(col, "col");
}

void TwoTetrahedral::CreateFiles()
{
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

    outfunc(v, "ori");

    outfunc(v, "ori2");

    matrix<double> col(nt + ns, 3);
    for (int i = 0; i < nt; i++)
    {
        col(i, 0) = 0.87;
        col(i, 1) = 0.94;
        col(i, 2) = 1.;
    }
    for (int i = nt; i < nt + ns; i++)
    {
        col(i, 0) = 1.;
        col(i, 1) = 0.5;
        col(i, 2) = 0.;
    }

    outfunc(col, "col");
}

#endif /* COMBOPATCHOUTPUT_CPP */
