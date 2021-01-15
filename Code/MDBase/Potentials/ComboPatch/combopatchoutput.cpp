#ifndef COMBOPATCHOUTPUT_CPP
#define COMBOPATCHOUTPUT_CPP

void ComboPatch::CreateFiles() {
    matrix<double> v(1, 3);
    v(0, 0) = 1.0;
    v(0, 1) = 0.0;
    v(0, 2) = 0.0;

    outfunc(v, "ori1");

    matrix<double> col(1, 3);
    col(0, 0) = 0.;
    col(0, 1) = 0.5;
    col(0, 2) = 0.5;

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
    col(0, 0) = 0.;
    col(0, 1) = 0.5;
    col(0, 2) = 0.5;

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

    outfunc(v, "ori1");

    matrix<double> v2(1,3);

    v2(0, 0) = 1.0;
    v2(0, 1) = 0.0;
    v2(0, 2) = 0.0;

    outfunc(v2, "ori2");

    matrix<double> col(nt+ns, 3);
    for(int i = 0 ; i < nt ; i++)
    {
    col(i, 0) = 0.0;
    col(i, 1) = 0.5;
    col(i, 2) = 0.5;
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

    outfunc(v, "ori1");

    outfunc(v, "ori2");

    matrix<double> col(nt + ns, 3);
    for (int i = 0; i < nt; i++)
    {
        col(i, 0) = 0.;
        col(i, 1) = 0.5;
        col(i, 2) = 0.5;
    }
    for (int i = nt; i < nt + ns; i++)
    {
        col(i, 0) = 1.;
        col(i, 1) = 0.5;
        col(i, 2) = 0.;
    }

    outfunc(col, "col");
}

void TwoTetrahedralAndSingle::CreateFiles()
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

    matrix<double> v2(1,3);

    v2(0, 0) = nx1;
    v2(0, 1) = ny1;
    v2(0, 2) = nz1;

    outfunc(v, "ori1");

    outfunc(v, "ori2");

    outfunc(v, "ori3");

    matrix<double> col(nf, 3);
    for (int i = 0; i < nt; i++)
    {
        col(i, 0) = 0.;
        col(i, 1) = 0.5;
        col(i, 2) = 0.5;
    }
    for (int i = nt; i < ns; i++)
    {
        col(i, 0) = 1.;
        col(i, 1) = 0.5;
        col(i, 2) = 0.;
    }
    for (int i = ns; i < nf; i++)
    {
        col(i, 0) = 1.;
        col(i, 1) = 0.5;
        col(i, 2) = 0.5;
    }

    outfunc(col, "col");
}

void TetrahedralWithBivalent::CreateFiles()
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

    outfunc(v, "ori1");

    matrix<double> v2(2, 3);

    v2(0, 0) = nx5;
    v2(0, 1) = ny5;
    v2(0, 2) = nz5;

    v2(1, 0) = nx6;
    v2(1, 1) = ny6;
    v2(1, 2) = nz6;

    outfunc(v2, "ori2");

    matrix<double> col(nt + nb, 3);
    for (int i = 0; i < nt; i++)
    {
        col(i, 0) = 0.0;
        col(i, 1) = 0.5;
        col(i, 2) = 0.5;
    }
    for (int i = nt; i < nt + nb; i++)
    {
        col(i, 0) = 1.;
        col(i, 1) = 0.5;
        col(i, 2) = 0.;
    }

    outfunc(col, "col");
}

#endif /* COMBOPATCHOUTPUT_CPP */
