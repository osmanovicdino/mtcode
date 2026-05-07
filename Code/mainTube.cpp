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
//  Particle design: 4 patches
//    Patch 0: north (+z pole)  — bonds to ring above
//    Patch 1: south (-z pole)  — bonds to ring below
//    Patch 2: lateral "left"   — bonds to counter-clockwise neighbour
//    Patch 3: lateral "right"  — bonds to clockwise neighbour
//
//  Patches 2 and 3 lie in the z=0 body-frame plane.  Their directions
//  are derived exactly from N_ring so that the tube closes on itself
//  for any ring size:
//
//    patch 2 body direction = (-sin(pi/N_ring),  cos(pi/N_ring), 0)
//    patch 3 body direction = (-sin(pi/N_ring), -cos(pi/N_ring), 0)
//
//  Derivation: particle i sits at lab azimuth phi_i = 2*pi*i/N_ring.
//  Its orientation matrix rotates body-x to the outward radial direction.
//  The unit bond vector from i to its CCW neighbour i+1 in the lab frame
//  is (-sin(phi_i + pi/N_ring), cos(phi_i + pi/N_ring), 0).  Applying
//  the body<-lab transform (multiply by R) cancels phi_i, leaving the
//  phi-independent body-frame directions above.  Patch 2 of particle i
//  therefore bonds perfectly with patch 3 of particle i+1, and vice
//  versa, for any N_ring.
//
//  Tube radius is set so the chord between ring neighbours equals 1.1
//  particle diameters (just above the WCA repulsive core):
//    R = 1.1 / (2 * sin(pi / N_ring))
//
//  Change N_ring to get wider or narrower tubes (6 = hexagonal,
//  10 = decagonal, etc.).
// ============================================================

int main(int argc, char **argv)
{
    srand(time(NULL));
    signal(SIGSEGV, handler);

    // ------------------------------------------------------------------
    // Geometry parameters  —  change N_ring to resize the tube
    // ------------------------------------------------------------------
    const int    N_ring   = 6;    // particles per ring (try 6, 8, 10, ...)
    const int    N_z      = 5;    // number of rings along z
    const int    N        = N_ring * N_z;
    const double chord    = 1.1;  // target neighbour distance (particle diameters)
    const double R        = chord / (2.0 * sin(pid / N_ring)); // tube radius
    const double dz       = chord; // axial ring spacing equals lateral spacing
    const double L        = 20.0; // cubic box side length

    // ------------------------------------------------------------------
    // Interaction parameters
    // ------------------------------------------------------------------
    const double axial_strength   = 15.0;
    const double lateral_strength = 15.0;
    const double range            = 1.2;  // interaction range (particle diameters)
    const double patchsize        = 0.5;  // max patch half-angle (radians)

    // ------------------------------------------------------------------
    // Simulation parameters
    // ------------------------------------------------------------------
    const int    runtime  = 10000000;
    const int    every    = 1000;
    const double kT       = 1.0;
    const double dt       = 0.005;
    const double viscosity= 1.0;
    const double hdradius = 0.5;

    // ==================================================================
    // Build the GeneralPatch potential
    // ==================================================================

    // One particle type, 4 patches, N particles total.
    vector1<int> vec1(1);
    vec1[0] = 4;

    vector1<int> numb(1);
    numb[0] = N;

    // Body-frame patch directions (4 rows x 3 cols).
    // Lateral patch directions encode the tube geometry via N_ring.
    const double lat_x = -sin(pid / N_ring);
    const double lat_y =  cos(pid / N_ring);

    matrix<double> patch_orient(4, 3);
    patch_orient(0, 0) = 0.0;    patch_orient(0, 1) = 0.0;     patch_orient(0, 2) =  1.0; // north
    patch_orient(1, 0) = 0.0;    patch_orient(1, 1) = 0.0;     patch_orient(1, 2) = -1.0; // south
    patch_orient(2, 0) = lat_x;  patch_orient(2, 1) =  lat_y;  patch_orient(2, 2) =  0.0; // left
    patch_orient(3, 0) = lat_x;  patch_orient(3, 1) = -lat_y;  patch_orient(3, 2) =  0.0; // right

    // Interaction matrix: 4x4 = 16 potentials, row index = patch_a * 4 + patch_b.
    // Only 4 bonding pairs are active:
    //   north(0) <-> south(1)  for axial bonding
    //   left(2)  <-> right(3)  for lateral bonding
    matrix<double> params(16, 3);
    for (int i = 0; i < 16; i++) {
        params(i, 0) = 0.0;
        params(i, 1) = range;
        params(i, 2) = patchsize;
    }
    params(0 * 4 + 1, 0) = axial_strength;   // north -> south
    params(1 * 4 + 0, 0) = axial_strength;   // south -> north
    params(2 * 4 + 3, 0) = lateral_strength; // left  -> right
    params(3 * 4 + 2, 0) = lateral_strength; // right -> left

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
