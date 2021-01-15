#ifndef SINGPATCH_CPP
#define SINGPATCH_CPP

SingPatch::SingPatch(double strr, double disss, double angg) : ComboPatch(1), ang(angg), dis(disss), str(strr)
{
    i1 = new int[2];
    i1[0] = 1;
    i1[1] = 0;
    p = &i1;

    mypot *pot1 = new mypot(nx, ny, nz, nx, ny, nz, str, dis, ang, 0.75);

    potential_bundle[0] = pot1->clone();

    delete pot1;
    //potential_bundle[0] = pot1;
    //vector1<potentialtheta3D *> pots(1);
}


int SingPatch::num_patches(const int &) { return 1; }

void SingPatch::UpdateIterator(const int &i, const int &j)
{
    //do nothing (only one patch)
}
void SingPatch::UpdateIteratorSafe(const int &i, const int &j, int **q)
{
    *q = i1;
}

int SingPatch::get_total_patches(const int &N) { return N; }
void SingPatch::which_patch(const int &i, const int &j, const int &potn, int &wpi, int &wpj)
{
    wpi = i;
    wpj = j;
}
void SingPatch::which_particle(const int &wpi, const int &wpj, int &i, int &j)
{
    i = wpi;
    j = wpj;
}
int SingPatch::which_potential(const int &i, const int &j, const int &wpi, const int &wpj)
{
    return 0;
}

void SingPatch::get_params(const int &i, const int &j, const int &potn, double &nxb1, double &nyb1, double &nzb1, double &nxb2, double &nyb2, double &nzb2, double &d12, double &ang12)
{
    nxb1 = nx;
    nyb1 = ny;
    nzb1 = nz;
    nxb2 = nx;
    nyb2 = ny;
    nzb2 = nz;
    d12 = dis;
    ang12 = ang;
}

#endif /* SINGPATCH_CPP */
