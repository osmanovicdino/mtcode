with these definitions of the torque found in potentialtheta3.h 

<pre><code>
            tix = (pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(2) * ny1 - un.gpcons(1) * nz1);
            tiy = (pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (-un.gpcons(2) * nx1 + un.gpcons(0) * nz1);
            tiz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(1) * nx1 - un.gpcons(0) * ny1);
            tjx = (pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(2) * ny2 + un.gpcons(1) * nz2);
            tjy = (pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (un.gpcons(2) * nx2 - un.gpcons(0) * nz2);
            tjz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(1) * nx2 + un.gpcons(0) * ny2);
 </code></pre>


it is found that this means that the attractive patches are no longer attractive to each other when in a particular configuration.

The configuration file for these torques to reproduce the bugs is given by the following main file:


<pre><code>

int main(int argc, char** argv) {
srand (time(NULL));


vector1<bool> pb(3,true);
cube geo(100.0,pb,3);

LangevinNVTR a(geo);


matrix<double> dat(2,3);

//for the given system we can rotate everything

matrix<double> Rot(3,3);

double thetax  = 0.0;

double thetay = pi/4.;

double thetaz = -pi/4.;

Rot(0, 0) = cos(thetay) * cos(thetaz);
Rot(0, 1) = -(cos(thetay) * sin(thetaz));
Rot(0, 2) = sin(thetay);
Rot(1, 0) = cos(thetaz) * sin(thetax) * sin(thetay) + cos(thetax) * sin(thetaz);
Rot(1, 1) = cos(thetax) * cos(thetaz) - sin(thetax) * sin(thetay) * sin(thetaz);
Rot(1, 2) = -(cos(thetay) * sin(thetax));
Rot(2, 0) = -(cos(thetax) * cos(thetaz) * sin(thetay)) + sin(thetax) * sin(thetaz);
Rot(2, 1) = cos(thetaz) * sin(thetax) + cos(thetax) * sin(thetay) * sin(thetaz);
Rot(2, 2) = cos(thetax) * cos(thetay);

// double nx = 1;
// double ny = 0;
// double nz = 0;
double theta = 0.2;
double phi = pi/2+0.2;

double nxb = cos(theta) * sin(phi);
double nyb = sin(theta) * sin(phi);
double nzb = cos(phi);

double nx = Rot(0, 0) * nxb + Rot(0, 1) * nyb + Rot(0,2) * nzb;

double ny = Rot(1, 0) * nxb + Rot(1, 1) * nyb + Rot(1, 2) * nzb;

double nz = Rot(2, 0) * nxb + Rot(2, 1) * nyb + Rot(2, 2) * nzb;

KernFrenkelOnePatch pot(nx,ny,nz,-nx,-ny,-nz, 100.0, 1.5, pi / 3., 0.75);

matrix<double> pos(2,3);
pos(0,0) = 50.0;
pos(0,1) = 50.0;
pos(0,2) = 50.0;
pos(1,0) = 52.0;
pos(1,1) = 50.0;
pos(1,2) = 50.0;

double pos1x = Rot(0, 0) * pos(0, 0) + Rot(0, 1) * pos(0, 1) + Rot(0, 2) * pos(0, 2);

double pos1y = Rot(1, 0) * pos(0, 0) + Rot(1, 1) * pos(0, 1) + Rot(1, 2) * pos(0, 2);

double pos1z = Rot(2, 0) * pos(0, 0) + Rot(2, 1) * pos(0, 1) + Rot(2, 2) * pos(0, 2);

double pos2x = Rot(0, 0) * pos(1, 0) + Rot(0, 1) * pos(1, 1) + Rot(0, 2) * pos(1, 2);

double pos2y = Rot(1, 0) * pos(1, 0) + Rot(1, 1) * pos(1, 1) + Rot(1, 2) * pos(1, 2);

double pos2z = Rot(2, 0) * pos(1, 0) + Rot(2, 1) * pos(1, 1) + Rot(2, 2) * pos(1, 2);

double tx = (pos1x - 50.0);

double ty = (pos1y - 50.0);

double tz = (pos1z - 50.0);

pos(0, 0) = pos1x - tx;
pos(0, 1) = pos1y - ty;
pos(0, 2) = pos1z - tz;
pos(1, 0) = pos2x - tx;
pos(1, 1) = pos2y - ty;
pos(1, 2) = pos2z - tz;

cout << Rot << endl;
cout << pos << endl;
pausel();

matrix<double> mom(2,3);
matrix<double> amo(2,3);

matrix<double> ori(2,9);

ori(0, 0) = 1;
ori(0, 1) = 0;
ori(0, 2) = 0;
ori(0, 3) = 0;
ori(0, 4) = 1;
ori(0, 5) = 0;
ori(0, 6) = 0;
ori(0, 7) = 0;
ori(0, 8) = 1;

ori(1, 0) = 1;
ori(1, 1) = 0;
ori(1, 2) = 0;
ori(1, 3) = 0;
ori(1, 4) = 1;
ori(1, 5) = 0;
ori(1, 6) = 0;
ori(1, 7) = 0;
ori(1, 8) = 1;


a.initialize(pos,mom,ori,amo);

matrix<double> F(2,3);
matrix<double> T(2, 3);
matrix<double> zeromatrix(2,3);
matrix<int> pairs(1,2);
pairs(0,0) = 0;
pairs(0,1) = 1;


WCAPotential wsa(1.0,1.0,0.0);


a.setdt(0.001);

double viscosity = 1.0;
double hdradius = 0.5;
double fac1 = 100.;
double fac2 = 1.;
a.setgamma(fac1*6*pi*viscosity*hdradius);
a.setgammar(fac2*8*pi*viscosity*hdradius*hdradius*hdradius);
a.setkT(0.0);

a.calculateforces(pairs, wsa);
a.calculate_forces_and_torques3D(pairs, pot, F, T);

matrix<double> tempT = T;

a.create_forces_and_torques_sphere(F, T);



ofstream myfile;
myfile.open("pos.csv");

ofstream myfile2;
myfile2.open("orientations.csv");

ofstream myfile3;
myfile3.open("patchpositions.csv");

ofstream myfile4;
myfile4.open("torques.csv");

ofstream myfile5;
myfile5.open("temperatures.csv");

for(int i = 0 ; i < 20000 ; i++ ) {

cout << i << endl;
a.measured_temperature(myfile5);

a.advancemom_halfstep(F,T);
a.advance_pos();
a.rotate();

F = a.calculateforces(pairs, wsa);
T = zeromatrix;
a.calculate_forces_and_torques3D(pairs, pot, F, T);

tempT = T;



a.create_forces_and_torques_sphere(F, T);


a.advancemom_halfstep(F,T);


// cout << a.getmom() << endl;
// cout << a.getangmom() << endl;
// pausel();
if(i%1 == 0) {

    matrix<double> orient = a.getorientation();
    matrix<double> pos = a.getdat();

    double qtemp0 = orient.gpcons(0, 0);
    double qtemp1 = orient.gpcons(0, 1);
    double qtemp2 = orient.gpcons(0, 2);
    double qtemp3 = orient.gpcons(0, 3);
    double qtemp4 = orient.gpcons(0, 4);
    double qtemp5 = orient.gpcons(0, 5);
    double qtemp6 = orient.gpcons(0, 6);
    double qtemp7 = orient.gpcons(0, 7);
    double qtemp8 = orient.gpcons(0, 8);

    double gtemp0 = orient.gpcons(1, 0);
    double gtemp1 = orient.gpcons(1, 1);
    double gtemp2 = orient.gpcons(1, 2);
    double gtemp3 = orient.gpcons(1, 3);
    double gtemp4 = orient.gpcons(1, 4);
    double gtemp5 = orient.gpcons(1, 5);
    double gtemp6 = orient.gpcons(1, 6);
    double gtemp7 = orient.gpcons(1, 7);
    double gtemp8 = orient.gpcons(1, 8);

    double nx1 = pot.nxb1 * qtemp0 + pot.nyb1 * qtemp3 + pot.nzb1 * qtemp6;
    double ny1 = pot.nxb1 * qtemp1 + pot.nyb1 * qtemp4 + pot.nzb1 * qtemp7;
    double nz1 = pot.nxb1 * qtemp2 + pot.nyb1 * qtemp5 + pot.nzb1 * qtemp8;

    double nx2 = pot.nxb2 * gtemp0 + pot.nyb2 * gtemp3 + pot.nzb2 * gtemp6;
    double ny2 = pot.nxb2 * gtemp1 + pot.nyb2 * gtemp4 + pot.nzb2 * gtemp7;
    double nz2 = pot.nxb2 * gtemp2 + pot.nyb2 * gtemp5 + pot.nzb2 * gtemp8;

    // cout << pos(0, 0) + 0.5 * nx1 << "," << pos(0, 1) + 0.5 * ny1 << "," << pos(0, 2) + 0.5 * nz1 << endl;
    // cout << pos(1, 0) + 0.5 * nx2 << "," << pos(1, 1) + 0.5 * ny2 << "," << pos(1, 2) + 0.5 * nz2 << endl;
    myfile3 << pos(0, 0) + 0.5 * nx1 << "," << pos(0, 1) + 0.5 * ny1 << "," << pos(0, 2) + 0.5 * nz1 << endl;
    myfile3 << pos(1, 0) + 0.5 * nx2 << "," << pos(1, 1) + 0.5 * ny2 << "," << pos(1, 2) + 0.5 * nz2 << endl;

    myfile <<= pos;
    myfile2 <<= orient;
    myfile4 <<= tempT;

}

}
myfile.close();
myfile2.close();
myfile3.close();
myfile4.close();



return 0;
}

</code></pre>

I will proceed by changing tx to -tx and see if that fixes the bug in all circumstances.

The following form seems to be working:

<pre><code>

            tix = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(2) * ny1 - un.gpcons(1) * nz1);
            tiy = (pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (-un.gpcons(2) * nx1 + un.gpcons(0) * nz1);
            tiz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(1) * nx1 - un.gpcons(0) * ny1);
            // cout << tx << " " << ty << " " << tz << endl;
            tjx = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(2) * ny2 + un.gpcons(1) * nz2);
            tjy = (pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (un.gpcons(2) * nx2 - un.gpcons(0) * nz2);
            tjz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(1) * nx2 + un.gpcons(0) * ny2);
      </code></pre>

Addendum: This also doesn't work.  The issue came from the interpretation of the propagation matrices in Sun et. al (2008). The matrix for R_y has a negative sign in the opposite direction. After updating genfullmat function to take account of this fact, the results now seem to behave for two particles. The torques will be given fully by:

<pre><code>

            tix = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(2) * ny1 - un.gpcons(1) * nz1);
            tiy = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (-un.gpcons(2) * nx1 + un.gpcons(0) * nz1);
            tiz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetai / (2 * thetam)) * d1 * cos(pi * thetaj / (2 * thetam)) * (un.gpcons(1) * nx1 - un.gpcons(0) * ny1);
            // cout << tx << " " << ty << " " << tz << endl;
            tjx = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(2) * ny2 + un.gpcons(1) * nz2);
            tjy = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (un.gpcons(2) * nx2 - un.gpcons(0) * nz2);
            tjz = -(pi / (2. * thetam)) * pot * 2 * f * sin(pi * thetaj / (2 * thetam)) * d2 * cos(pi * thetai / (2 * thetam)) * (-un.gpcons(1) * nx2 + un.gpcons(0) * ny2);
      </code></pre>      