#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>
#include <random>
#include <algorithm>
#include <parallel/algorithm>
#include <unistd.h>
#include <execinfo.h>
#include <signal.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

using namespace std;

#include "Basic/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"
#include "MDBase/MD.h"
#include "MDBase/Langevin.h"
#include "MDBase/Rotational/LangevinR.h"

// ============================================================
//  Hollow Tube Simulation
//
//  Particle design: 5 patches
//    Patch 0: north (+z pole)
//    Patch 1: south (-z pole)
//    Patches 2,3,4: three lateral patches at polar angle theta from z,
//                   separated 120 degrees apart in azimuth.
//
//  For a ring of N_ring = 6 particles (hexagonal cross-section) with
//  equatorial lateral patches (theta = pi/2), the geometry works out
//  exactly: patch 3 of particle i bonds to patch 4 of its clockwise
//  neighbour, and vice versa.  The axial patches bond particles in
//  adjacent rings directly above / below.
//
//  Initial positions place particles on a cylinder:
//    (x, y, z) = (R cos phi_i + L/2,  R sin phi_i + L/2,  j * dz + z0)
//  with phi_i = 2 pi i / N_ring.  Each particle's orientation matrix is
//  set so that its body-frame x-axis points radially outward, keeping
//  the lateral patches aligned with their intended bond partners from
//  the very first step.
// ============================================================

int main(int argc, char **argv)
{
    srand(time(NULL));
    signal(SIGSEGV, handler);

    // ------------------------------------------------------------------
    // Geometry parameters
    // ------------------------------------------------------------------
    const int    N_ring   = 6;           // particles per ring (hexagonal)
    const int    N_z      = 5;           // number of rings along z
    const int    N        = N_ring * N_z;// total particles
    const double R        = 1.5;         // tube radius (particle diameters)
    const double dz       = 1.2;         // z-spacing between rings
    const double theta    = pid / 2.0;   // polar angle from z for lateral patches
    const double L        = 20.0;        // cubic box side length

    // ------------------------------------------------------------------
    // Interaction parameters
    // ------------------------------------------------------------------
    const double axial_strength   = 15.0; // north-south bond energy
    const double lateral_strength = 15.0; // lateral bond energy
    const double range            = 1.2;  // interaction range (particle diameters)
    const double patchsize        = 0.6;  // max patch half-angle (radians, ~34 deg)

    // ------------------------------------------------------------------
    // Simulation parameters
    // ------------------------------------------------------------------
    const int    runtime  = 10000000;
    const int    every    = 1000;
    const double kT       = 1.0;
    const double dt       = 0.005;
    const double viscosity= 1.0;
    const double hdradius = 0.5;         // hydrodynamic radius

    // ==================================================================
    // Build the GeneralPatch potential
    // ==================================================================

    // One particle type with 5 patches; N particles total.
    vector1<int> vec1(1);
    vec1[0] = 5;

    vector1<int> numb(1);
    numb[0] = N;

    // Patch body-frame directions (5 rows x 3 cols)
    // Row 0: north pole  (0, 0,  1)
    // Row 1: south pole  (0, 0, -1)
    // Rows 2-4: lateral patches at polar angle theta, azimuth 0, 120, 240 deg
    matrix<double> patch_orient(5, 3);

    patch_orient(0, 0) = 0.0;   patch_orient(0, 1) = 0.0;                   patch_orient(0, 2) =  1.0;
    patch_orient(1, 0) = 0.0;   patch_orient(1, 1) = 0.0;                   patch_orient(1, 2) = -1.0;

    patch_orient(2, 0) = sin(theta);
    patch_orient(2, 1) = 0.0;
    patch_orient(2, 2) = cos(theta);

    patch_orient(3, 0) = sin(theta) * cos(2.0 * pid / 3.0);
    patch_orient(3, 1) = sin(theta) * sin(2.0 * pid / 3.0);
    patch_orient(3, 2) = cos(theta);

    patch_orient(4, 0) = sin(theta) * cos(4.0 * pid / 3.0);
    patch_orient(4, 1) = sin(theta) * sin(4.0 * pid / 3.0);
    patch_orient(4, 2) = cos(theta);

    // Interaction parameters: 5x5 = 25 potentials (cols: strength, range, patchsize)
    // Indexing: potential for patch a of particle 1 with patch b of particle 2 is row a*5+b.
    matrix<double> params(25, 3);
    for (int i = 0; i < 25; i++) {
        params(i, 0) = 0.0;      // default: no interaction
        params(i, 1) = range;
        params(i, 2) = patchsize;
    }

    // Axial bonding: north (patch 0) with south (patch 1) and vice versa
    params(0 * 5 + 1, 0) = axial_strength;  // north -> south
    params(1 * 5 + 0, 0) = axial_strength;  // south -> north

    // Lateral bonding: patches 2, 3, 4 interact with each other
    // For N_ring=6, the active bonding pairs are 3<->4 and 4<->3;
    // keeping all lateral cross-terms on lets the structure self-heal.
    for (int a = 2; a < 5; a++) {
        for (int b = 2; b < 5; b++) {
            params(a * 5 + b, 0) = lateral_strength;
        }
    }

    GeneralPatch c(vec1, numb, params, patch_orient);

    // ==================================================================
    // Initial tube positions and orientations
    // ==================================================================

    // Positions: particle k = ring j * N_ring + i
    //   sits at (R cos phi_i, R sin phi_i, j*dz) centred in the box.
    // Orientations: rotation matrix R such that n_lab = R^T * n_body,
    //   mapping body-x to the outward radial direction at azimuth phi_i.
    //     R = [ cos(phi)   sin(phi)   0 ]
    //         [-sin(phi)   cos(phi)   0 ]
    //         [   0          0        1 ]
    // Stored row-major: orient(k, 0..8) = R[0][0..2], R[1][0..2], R[2][0..2].

    matrix<double> pos(N, 3);
    matrix<double> moms(N, 3);       // zero initial momenta
    matrix<double> orients(N, 9);
    matrix<double> angmoms(N, 3);    // zero initial angular momenta

    double z0 = L / 2.0 - (N_z - 1) * dz / 2.0;  // centre tube in box along z

    for (int j = 0; j < N_z; j++) {
        for (int i = 0; i < N_ring; i++) {
            int k = j * N_ring + i;
            double phi = 2.0 * pid * i / N_ring;

            pos(k, 0) = L / 2.0 + R * cos(phi);
            pos(k, 1) = L / 2.0 + R * sin(phi);
            pos(k, 2) = z0 + j * dz;

            // Orientation matrix (row-major 3x3)
            orients(k, 0) =  cos(phi);  orients(k, 1) = sin(phi);  orients(k, 2) = 0.0;
            orients(k, 3) = -sin(phi);  orients(k, 4) = cos(phi);  orients(k, 5) = 0.0;
            orients(k, 6) = 0.0;        orients(k, 7) = 0.0;        orients(k, 8) = 1.0;
        }
    }

    // ==================================================================
    // Langevin integrator setup
    // ==================================================================

    vector1<bool> pb(3, true);  // periodic boundaries in all directions
    cube geo(L, pb, 3);

    LangevinNVTR obj(geo);
    obj.initialize(pos, moms, orients, angmoms);

    obj.setdt(dt);
    obj.setgamma(6.0 * pid * viscosity * hdradius);
    obj.setgammar(8.0 * pid * viscosity * hdradius * hdradius * hdradius);
    obj.setkT(kT);

    // ==================================================================
    // Run loop
    // ==================================================================

    int num_cells = (int)floor(L / 4.0);
    int ccc;
    matrix<int> boxes = geo.generate_boxes_relationships(num_cells, ccc);
    matrix<int> *pairs = obj.calculatepairs(boxes, 3.5);

    // Soft repulsive core to prevent overlaps
    WCAPotential wsa(3.0, 1.0, 0.0);

    matrix<double> F(N, 3);
    matrix<double> T(N, 3);
    matrix<double> RT(N, 6);

    // First force/torque evaluation
    F = obj.calculateforces(*pairs, wsa);
    obj.calculate_forces_and_torques3D(*pairs, c, F, T);
    generate_uniform_random_matrix(RT);
    obj.create_forces_and_torques_sphere(F, T, RT);

    int number_of_digits = 0;
    {
        int tf = (runtime / every) + 1;
        do { ++number_of_digits; tf /= 10; } while (tf);
    }

    for (int i = 0; i < runtime; i++) {
        if (i > 0 && i % 20 == 0) {
            delete pairs;
            pairs = obj.calculatepairs(boxes, 3.5);
        }

        obj.advancemom_halfstep(F, T);
        obj.advance_pos();
        obj.rotate();

        F = obj.calculateforces(*pairs, wsa);
        T.reset(0.0);
        obj.calculate_forces_and_torques3D(*pairs, c, F, T);
        generate_uniform_random_matrix(RT);
        obj.create_forces_and_torques_sphere(F, T, RT);
        obj.advancemom_halfstep(F, T);

        if (i % every == 0) {
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << (i / every);
            outfunc(obj.getdat(),         "tube_pos_i=" + ss.str() + ".csv");
            outfunc(obj.getorientation(), "tube_ori_i=" + ss.str() + ".csv");
        }
    }

    delete pairs;
    return 0;
}
