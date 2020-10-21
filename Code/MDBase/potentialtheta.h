#ifndef POTENTIALTHETA_H
#define POTENTIALTHETA_H


#include "../Basic/basic.h"

struct potentialtheta3D
{            //this is the base constructor for all potentials.
    bool dl; //does this interaction decay at large distances?
    double interaction_distance;

    virtual potentialtheta3D *clone() const = 0;

    // virtual double energy(double, double) = 0;
    // virtual double force(vector1<double> &pos, const matrix<double> &A1,const matrix<double> &A2 ) = 0;
    virtual void force_and_torque(const vector1<double> &un, double rij, const matrix<double> &orient, int i, int j, double &, double &, double &, double &, double &, double &, double &, double &, double &) = 0;
    // virtual double torque(double,double) = 0;
    //virtual double force2mdx(double, double) = 0; //the factor as a function of distance*squared which when multiplied with dx will give the total force

    virtual void setparameters(const vector1<double> &) = 0;
    // virtual void printparameters() = 0;
    // virtual void printparameters(ofstream&) = 0;
};

struct KernFrenkelOnePatch : potentialtheta3D {
    // Q lab -> body
    // Q^T body -> lab

    double nxb1; //positon to the patch relative to centre of mass of the sphere
    double nyb1;
    double nzb1;
 
    double nxb2;
    double nyb2;
    double nzb2;

    double att;
    double dis;
    double thetam;

    double v;


    KernFrenkelOnePatch(double nxx1, double nyy1, double nzz1, double nxx2, double nyy2, double nzz2, double attt, double diss, double thetamm, double vv)
    {

        nxb1 = nxx1;
        nyb1 = nyy1;
        nzb1 = nzz1;
        nxb2 = nxx2;
        nyb2 = nyy2;
        nzb2 = nzz2;
        if (abs(SQR(nxb1) + SQR(nyb2) + SQR(nzb2) - 1) > 1E-5)
            error("patch potential should be called with unit vector (vector not normalized)");
        if (abs(SQR(nxb2) + SQR(nyb2) + SQR(nzb2) - 1) > 1E-5)
            error("patch potential should be called with unit vector (vector not normalized)");

        dl = false;
        dis = diss;
        att = attt;
        interaction_distance =  2*dis;
        thetam = thetamm;
        //v = vv;
        v =  vv;
    }

    void force_and_torque(const vector1<double> &un, double rij, const matrix<double> &orient, int i, int j, double &fx, double &fy, double &fz, double &tix, double &tiy, double &tiz, double &tjx, double &tjy, double &tjz)
    {
        //dpos is the vector of differences of distance in each dimension
        if(rij > 2.5*(dis+v) ) {
            fx = 0.0;
            fy = 0.0;
            fz = 0.0;
            tix = 0.0;
            tiy = 0.0;
            tiz = 0.0;
            tjx = 0.0;
            tjy = 0.0;
            tjz = 0.0;
        }
        else{
            double qtemp0 = orient.gpcons(i, 0);
            double qtemp1 = orient.gpcons(i, 1);
            double qtemp2 = orient.gpcons(i, 2);
            double qtemp3 = orient.gpcons(i, 3);
            double qtemp4 = orient.gpcons(i, 4);
            double qtemp5 = orient.gpcons(i, 5);
            double qtemp6 = orient.gpcons(i, 6);
            double qtemp7 = orient.gpcons(i, 7);
            double qtemp8 = orient.gpcons(i, 8);

            double gtemp0 = orient.gpcons(j, 0);
            double gtemp1 = orient.gpcons(j, 1);
            double gtemp2 = orient.gpcons(j, 2);
            double gtemp3 = orient.gpcons(j, 3);
            double gtemp4 = orient.gpcons(j, 4);
            double gtemp5 = orient.gpcons(j, 5);
            double gtemp6 = orient.gpcons(j, 6);
            double gtemp7 = orient.gpcons(j, 7);
            double gtemp8 = orient.gpcons(j, 8);

            double nx1 = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
            double ny1 = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
            double nz1 = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;

            double nx2 = nxb2 * gtemp0 + nyb2 * gtemp3 + nzb2 * gtemp6;
            double ny2 = nxb2 * gtemp1 + nyb2 * gtemp4 + nzb2 * gtemp7;
            double nz2 = nxb2 * gtemp2 + nyb2 * gtemp5 + nzb2 * gtemp8;

            double argthetai = -(nx1 * un.gpcons(0) + ny1 * un.gpcons(1) + nz1 * un.gpcons(2));
            double argthetaj = (nx2 * un.gpcons(0) + ny2 * un.gpcons(1) + nz2 * un.gpcons(2));

            double f;




            if(argthetai >cos(thetam) && argthetaj > cos(thetam)) {

            double thetai = acos(argthetai);

            double thetaj = acos(argthetaj);

            f = cos(pi * thetai / (2. * thetam)) * cos(pi * thetaj / (2. * thetam));

            double fac = (dis/(rij+v*dis));
            double fac2 = SQR(fac);
            double fac6 = CUB(fac2);
            double fac12 = SQR(fac6);

            double potf = ((24*att)/dis)*(2*fac*fac12-fac*fac6);

            double pot =  4*att*((fac12)-(fac6));




            fx =  SQR(f) * potf * un.gpcons(0);
            fy =  SQR(f) * potf * un.gpcons(1);
            fz =  SQR(f) * potf * un.gpcons(2);

            //cout << fx << " " << fy << " " << fz << endl;

            double d1;
            if(abs(SQR(argthetai)-1)<1E-5) {
                d1 = 0.0;
            }
            else{
                d1 = -1. / sqrt(1 - SQR(argthetai));
            }
            double d2;
            if (abs(SQR(argthetaj) - 1) < 1E-5)
            {
                d2 = 0.0;
            }
            else {
                d2 = -1. / sqrt(1 - SQR(argthetaj));
            }

            // double factor = ((rij/dis)-SQR(rij/dis));
           // cout << "DONE" << endl;

            double temp1  = -pot*2*f*((pi/(2*thetam))*sin(pi * thetai / (2. * thetam))*d1*cos((pi*thetaj)/(2*thetam)));

            double temp2  = -pot*2*f*((pi/(2*thetam))*sin(pi * thetaj / (2. * thetam))*d2*cos((pi*thetai)/(2*thetam)));

           // cout << fx << " " << fy << " " << fz << endl;

            fx += (temp1/rij)*(nx1 + (argthetai) * un.gpcons(0));
            fy += (temp1/rij)*(ny1 + (argthetai) * un.gpcons(1));
            fz += (temp1/rij)*(nz1 + (argthetai) * un.gpcons(2));
            //cout << temp1 * (-nx1 / rij + (argthetai / SQR(rij)) * un.gpcons(0)) << " " << temp1 * (-ny1 / rij + (argthetai / SQR(rij)) * un.gpcons(1)) << " " << fz << endl;

            // cout << fx << " " << fy << " " << fz << endl;
            fx += -(temp2/rij)*(nx2  - (argthetaj) * un.gpcons(0));
            fy += -(temp2/rij)*(ny2  - (argthetaj) * un.gpcons(1));
            fz += -(temp2/rij)*(nz2  - (argthetaj) * un.gpcons(2));

            // cout << fx << " " << fy << " " << fz << endl;

            tix = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(2) * ny1 - un.gpcons(1) * nz1);
            tiy = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (-un.gpcons(2) * nx1 + un.gpcons(0) * nz1);
            tiz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(1) * nx1 - un.gpcons(0) * ny1);
            // cout << tx << " " << ty << " " << tz << endl;
            tjx = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(2) * ny2 + un.gpcons(1) * nz2);
            tjy = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (un.gpcons(2) * nx2 - un.gpcons(0) * nz2);
            tjz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(1) * nx2 + un.gpcons(0) * ny2);
            // cout << tx << " " << ty << " " << tz << endl;
            // cout << ttx << " " << tty << " " << ttz << endl;
            // //cout << fx << " " << fy << " " << fz << endl;

        //     double rf1x = -rij * un.gpcons(2) * fy + rij * un.gpcons(1) * fz;
        //     double rf1y =  rij * un.gpcons(2) * fx - rij * un.gpcons(0) * fz;
        //     double rf1z =  -rij * un.gpcons(1)*fx + rij * un.gpcons(0)*fy;


        //     cout << rf1x << " " << rf1y << " " << rf1z << endl;
        // cout << fx << " " << fy <<  " " << fz << endl;
        //     cout << rf1x + tix + tjx << "  " << rf1y + tiy + tjy << "  " << rf1z + tiz + tjz << endl;
        //     pausel();
            
            // pausel();

            //if (abs(tx - ttx) > 1E-10) error("not the same");
            // pausel();
            // cout << fx << " " << fy << " " << fz << endl;
            // pausel();
            
            // cout << un.gpcons(0) << " " << un.gpcons(1) << " " << un.gpcons(2) << endl;

            // cout << nx1 << " " << ny1 << " " << nz1 << endl;
            // cout << "all elements" << endl;
            // cout << temp1 << " " << temp2 << endl;
            // cout << -(pi * att * dis / (4 * thetam)) << endl;
            // cout << factor << endl;
            // cout << 2 * f << endl; 
            // cout << sin(pi * thetai / (2 * thetam)) << endl;
            // cout <<  d1 << endl;  
            // cout << cos(pi * thetaj / (2 * thetam)) << endl;
            // cout << (un.gpcons(2) * ny1 - un.gpcons(1) * nz1) << " " << (-un.gpcons(2) * nx1 + un.gpcons(0) * nz1) << " " << (un.gpcons(1) * nx1 - un.gpcons(0) * ny1) << endl;
            // cout << tx << " " << ty << " " << tz << endl;
            // pausel();

            // tx = (pi * att * dis / (4 * thetam)) * factor * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(0));
            // ty = (pi * att * dis / (4 * thetam)) * factor * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(1));
            // tz = (pi * att * dis / (4 * thetam)) * factor * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(2));
            }
            else {
                f = 0.;
                fx = 0.0;
                fy = 0.0;
                fz = 0.0;
                tix = 0.0;
                tiy = 0.0;
                tiz = 0.0;
                tjx = 0.0;
                tjy = 0.0;
                tjz = 0.0;
            }
        }
    }

    void setparameters(const vector1<double> &param) {
        nxb1 = param.gpcons(0);
        nyb1 = param.gpcons(1);
        nzb1 = param.gpcons(2);
        nxb2 = param.gpcons(3);
        nyb2 = param.gpcons(4);
        nzb2 = param.gpcons(5);
        if (abs(SQR(nxb1) + SQR(nyb2) + SQR(nzb2) - 1) > 1E-5)
            error("patch potential should be called with unit vector (vector not normalized)");
        if (abs(SQR(nxb2) + SQR(nyb2) + SQR(nzb2) - 1) > 1E-5)
            error("patch potential should be called with unit vector (vector not normalized)");

        dis = param.gpcons(6);
        att = param.gpcons(7);
        thetam = param.gpcons(8);
        v = param.gpcons(9);
    }

    KernFrenkelOnePatch *clone() const
    {
        return new KernFrenkelOnePatch(*this);
    }


};

struct KernFrenkelOnePatch2 : potentialtheta3D
{
    // Q lab -> body
    // Q^T body -> lab

    double nxb1; //positon to the patch relative to centre of mass of the sphere
    double nyb1;
    double nzb1;

    double nxb2;
    double nyb2;
    double nzb2;

    double att;
    double dis;
    double thetam;

    double v;

    KernFrenkelOnePatch2(double nxx1, double nyy1, double nzz1, double nxx2, double nyy2, double nzz2, double attt, double diss, double thetamm, double vv)
    {

        nxb1 = nxx1;
        nyb1 = nyy1;
        nzb1 = nzz1;
        nxb2 = nxx2;
        nyb2 = nyy2;
        nzb2 = nzz2;
        if (abs(SQR(nxb1) + SQR(nyb1) + SQR(nzb1) - 1) > 1E-5) {
            cout << SQR(nxb1) + SQR(nyb1) + SQR(nzb1) << endl;
            error("patch potential should be called with unit vector (vector not normalized)");
        }
        if (abs(SQR(nxb2) + SQR(nyb2) + SQR(nzb2) - 1) > 1E-5) {
            cout << SQR(nxb2) + SQR(nyb2) + SQR(nzb2) << endl;
            error("patch potential should be called with unit vector (vector not normalized)");
        }
        dl = false;
        dis = diss;
        att = attt;
        interaction_distance = 2 * dis;
        thetam = thetamm;
        //v = vv;
        v = vv;
    }

    void force_and_torque(const vector1<double> &un, double rij, const matrix<double> &orient, int i, int j, double &fx, double &fy, double &fz, double &tix, double &tiy, double &tiz, double &tjx, double &tjy, double &tjz)
    {
        //dpos is the vector of differences of distance in each dimension
        if (rij > 2.5 * (dis + v))
        {
            fx = 0.0;
            fy = 0.0;
            fz = 0.0;
            tix = 0.0;
            tiy = 0.0;
            tiz = 0.0;
            tjx = 0.0;
            tjy = 0.0;
            tjz = 0.0;
        }
        else
        {
            double qtemp0 = orient.gpcons(i, 0);
            double qtemp1 = orient.gpcons(i, 1);
            double qtemp2 = orient.gpcons(i, 2);
            double qtemp3 = orient.gpcons(i, 3);
            double qtemp4 = orient.gpcons(i, 4);
            double qtemp5 = orient.gpcons(i, 5);
            double qtemp6 = orient.gpcons(i, 6);
            double qtemp7 = orient.gpcons(i, 7);
            double qtemp8 = orient.gpcons(i, 8);

            double gtemp0 = orient.gpcons(j, 0);
            double gtemp1 = orient.gpcons(j, 1);
            double gtemp2 = orient.gpcons(j, 2);
            double gtemp3 = orient.gpcons(j, 3);
            double gtemp4 = orient.gpcons(j, 4);
            double gtemp5 = orient.gpcons(j, 5);
            double gtemp6 = orient.gpcons(j, 6);
            double gtemp7 = orient.gpcons(j, 7);
            double gtemp8 = orient.gpcons(j, 8);

            double nx1 = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
            double ny1 = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
            double nz1 = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;

            double nx2 = nxb2 * gtemp0 + nyb2 * gtemp3 + nzb2 * gtemp6;
            double ny2 = nxb2 * gtemp1 + nyb2 * gtemp4 + nzb2 * gtemp7;
            double nz2 = nxb2 * gtemp2 + nyb2 * gtemp5 + nzb2 * gtemp8;

            double argthetai = -(nx1 * un.gpcons(0) + ny1 * un.gpcons(1) + nz1 * un.gpcons(2));
            double argthetaj = (nx2 * un.gpcons(0) + ny2 * un.gpcons(1) + nz2 * un.gpcons(2));

            double f;

            if (argthetai > cos(thetam) && argthetaj > cos(thetam))
            {

                double thetai = acos(argthetai);

                double thetaj = acos(argthetaj);

                f = cos(pi * thetai / (2. * thetam)) * cos(pi * thetaj / (2. * thetam));

                double fac = (rij/dis);
                double fac2 = SQR(fac);
                double fac4 = SQR(fac2);
                double fac8 = SQR(fac4);
                double fac10 = fac2*fac8;

                double expf = exp(-0.5 * fac10);


                double pot = - att * expf;

                double potf = -5. * att * fac8 * (fac / dis) * expf;
                //double potf = ((24 * att) / dis) * (2 * fac * fac12 - fac * fac6);

                

                fx = SQR(f) * potf * un.gpcons(0);
                fy = SQR(f) * potf * un.gpcons(1);
                fz = SQR(f) * potf * un.gpcons(2);

                //cout << fx << " " << fy << " " << fz << endl;

                double d1;
                if (abs(SQR(argthetai) - 1) < 1E-5)
                {
                    d1 = 0.0;
                }
                else
                {
                    d1 = -1. / sqrt(1 - SQR(argthetai));
                }
                double d2;
                if (abs(SQR(argthetaj) - 1) < 1E-5)
                {
                    d2 = 0.0;
                }
                else
                {
                    d2 = -1. / sqrt(1 - SQR(argthetaj));
                }

                // double factor = ((rij/dis)-SQR(rij/dis));
                // cout << "DONE" << endl;

                double temp1 = -pot * 2 * f * ((pi / (2 * thetam)) * sin(pi * thetai / (2. * thetam)) * d1 * cos((pi * thetaj) / (2 * thetam)));

                double temp2 = -pot * 2 * f * ((pi / (2 * thetam)) * sin(pi * thetaj / (2. * thetam)) * d2 * cos((pi * thetai) / (2 * thetam)));

                // cout << fx << " " << fy << " " << fz << endl;

                fx += (temp1 / rij) * (nx1 + (argthetai)*un.gpcons(0));
                fy += (temp1 / rij) * (ny1 + (argthetai)*un.gpcons(1));
                fz += (temp1 / rij) * (nz1 + (argthetai)*un.gpcons(2));
                //cout << temp1 * (-nx1 / rij + (argthetai / SQR(rij)) * un.gpcons(0)) << " " << temp1 * (-ny1 / rij + (argthetai / SQR(rij)) * un.gpcons(1)) << " " << fz << endl;

                // cout << fx << " " << fy << " " << fz << endl;
                fx += -(temp2 / rij) * (nx2 - (argthetaj)*un.gpcons(0));
                fy += -(temp2 / rij) * (ny2 - (argthetaj)*un.gpcons(1));
                fz += -(temp2 / rij) * (nz2 - (argthetaj)*un.gpcons(2));

                // cout << fx << " " << fy << " " << fz << endl;

                tix = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(2) * ny1 - un.gpcons(1) * nz1);
                tiy = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (-un.gpcons(2) * nx1 + un.gpcons(0) * nz1);
                tiz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(1) * nx1 - un.gpcons(0) * ny1);
                // cout << tx << " " << ty << " " << tz << endl;
                tjx = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(2) * ny2 + un.gpcons(1) * nz2);
                tjy = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (un.gpcons(2) * nx2 - un.gpcons(0) * nz2);
                tjz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(1) * nx2 + un.gpcons(0) * ny2);
                // cout << tx << " " << ty << " " << tz << endl;
                // cout << ttx << " " << tty << " " << ttz << endl;
                // //cout << fx << " " << fy << " " << fz << endl;

                //     double rf1x = -rij * un.gpcons(2) * fy + rij * un.gpcons(1) * fz;
                //     double rf1y =  rij * un.gpcons(2) * fx - rij * un.gpcons(0) * fz;
                //     double rf1z =  -rij * un.gpcons(1)*fx + rij * un.gpcons(0)*fy;

                //     cout << rf1x << " " << rf1y << " " << rf1z << endl;
                // cout << fx << " " << fy <<  " " << fz << endl;
                //     cout << rf1x + tix + tjx << "  " << rf1y + tiy + tjy << "  " << rf1z + tiz + tjz << endl;
                //     pausel();

                // pausel();

                //if (abs(tx - ttx) > 1E-10) error("not the same");
                // pausel();
                // cout << fx << " " << fy << " " << fz << endl;
                // pausel();

                // cout << un.gpcons(0) << " " << un.gpcons(1) << " " << un.gpcons(2) << endl;

                // cout << nx1 << " " << ny1 << " " << nz1 << endl;
                // cout << "all elements" << endl;
                // cout << temp1 << " " << temp2 << endl;
                // cout << -(pi * att * dis / (4 * thetam)) << endl;
                // cout << factor << endl;
                // cout << 2 * f << endl;
                // cout << sin(pi * thetai / (2 * thetam)) << endl;
                // cout <<  d1 << endl;
                // cout << cos(pi * thetaj / (2 * thetam)) << endl;
                // cout << (un.gpcons(2) * ny1 - un.gpcons(1) * nz1) << " " << (-un.gpcons(2) * nx1 + un.gpcons(0) * nz1) << " " << (un.gpcons(1) * nx1 - un.gpcons(0) * ny1) << endl;
                // cout << tx << " " << ty << " " << tz << endl;
                // pausel();

                // tx = (pi * att * dis / (4 * thetam)) * factor * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(0));
                // ty = (pi * att * dis / (4 * thetam)) * factor * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(1));
                // tz = (pi * att * dis / (4 * thetam)) * factor * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(2));
            }
            else
            {
                f = 0.;
                fx = 0.0;
                fy = 0.0;
                fz = 0.0;
                tix = 0.0;
                tiy = 0.0;
                tiz = 0.0;
                tjx = 0.0;
                tjy = 0.0;
                tjz = 0.0;
            }
        }
    }

    void setparameters(const vector1<double> &param)
    {
        nxb1 = param.gpcons(0);
        nyb1 = param.gpcons(1);
        nzb1 = param.gpcons(2);
        nxb2 = param.gpcons(3);
        nyb2 = param.gpcons(4);
        nzb2 = param.gpcons(5);
        if (abs(SQR(nxb1) + SQR(nyb2) + SQR(nzb2) - 1) > 1E-5)
            error("patch potential should be called with unit vector (vector not normalized)");
        if (abs(SQR(nxb2) + SQR(nyb2) + SQR(nzb2) - 1) > 1E-5)
            error("patch potential should be called with unit vector (vector not normalized)");

        dis = param.gpcons(6);
        att = param.gpcons(7);
        thetam = param.gpcons(8);
        v = param.gpcons(9);
    }

    KernFrenkelOnePatch2 *clone() const
    {
        return new KernFrenkelOnePatch2(*this);
    }
};

#endif