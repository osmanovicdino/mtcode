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

// ==========================================================================
//  Two hollow tubes + binding-site lattice + star polymers
// ==========================================================================
//
//  TUBE PARTICLES (type 0) — 5 patches each
//    patch 0:  north (+z body axis)            -> axial bond to ring above
//    patch 1:  south (-z body axis)            -> axial bond to ring below
//    patch 2:  lateral "left"                  -> bond to CCW neighbour
//    patch 3:  lateral "right"                 -> bond to CW neighbour
//    patch 4:  outward radial (body +x)        -> BINDING SITE (lattice)
//
//  Two tubes of the same type are placed side by side along x.  Each has
//  N_ring particles per ring and N_z rings, so each tube exposes a
//  cylindrical lattice of N_ring x N_z binding sites with lateral spacing
//  ~chord and axial spacing dz.
//
//  STAR POLYMERS (type 1) — M_arms patches each
//    Arm-end patches are distributed approximately uniformly on a sphere
//    using a Fibonacci spiral.  Each arm-end can bond to a tube binding
//    site (patch 4 of type 0).
//
//    LIMITATION: this is a *rigid* M-patch sphere, not a true bead-spring
//    star polymer.  All M arm ends rotate as one body; flexible arms
//    would require bead-spring chains (possible in this codebase but a
//    much larger change).  This rigid model captures multivalent binding
//    geometry but loses the entropic flexibility of real polymer arms.
//
//  INTERACTIONS (only the four pairs below are non-zero):
//    tube  patch 0  <->  tube  patch 1    (axial,   strength axial_strength)
//    tube  patch 2  <->  tube  patch 3    (lateral, strength lateral_strength)
//    tube  patch 4  <->  star  patch *    (binding, strength binding_strength)
//    everything else: zero
// ==========================================================================

int main(int argc, char **argv)
{
    srand(time(NULL));
    signal(SIGSEGV, handler);

    // ----------------------------------------------------------------------
    // Tube geometry  (change N_ring to resize, N_z for length, N_tubes = 2)
    // ----------------------------------------------------------------------
    const int    N_ring    = 6;      // particles per ring (also = binding sites per slice)
    const int    N_z       = 6;      // number of rings per tube
    const int    N_tubes   = 2;
    const int    N_tube_each = N_ring * N_z;
    const int    N_tube    = N_tubes * N_tube_each;
    const double chord     = 1.1;    // lateral neighbour distance (particle diameters)
    const double R_tube    = chord / (2.0 * sin(pid / N_ring));
    const double dz        = chord;  // axial ring spacing = "distance N to another ring"

    // ----------------------------------------------------------------------
    // Star polymer parameters
    // ----------------------------------------------------------------------
    const int    N_star    = 30;     // number of free star polymers
    const int    M_arms    = 6;      // arms (= arm-end patches) per star

    // ----------------------------------------------------------------------
    // Box and bulk parameters
    // ----------------------------------------------------------------------
    const double tube_sep  = 3.0 * R_tube + 4.0; // distance between the two tube axes
    const double L         = 30.0;               // cubic box side length

    // ----------------------------------------------------------------------
    // Interaction strengths
    // ----------------------------------------------------------------------
    const double axial_strength   = 15.0;
    const double lateral_strength = 15.0;
    const double binding_strength = 8.0;
    const double range            = 1.2;   // interaction range (particle diameters)
    const double patchsize        = 0.5;   // max patch half-angle (radians)

    // ----------------------------------------------------------------------
    // Simulation parameters
    // ----------------------------------------------------------------------
    const int    runtime  = 10000000;
    const int    every    = 1000;
    const double kT       = 1.0;
    const double dt       = 0.005;
    const double viscosity= 1.0;
    const double hdradius = 0.5;

    const int N = N_tube + N_star;
    cout << "Setup: " << N_tubes << " tubes x " << N_tube_each
         << " tube particles + " << N_star << " star polymers (M=" << M_arms
         << ") = " << N << " particles total" << endl;
    cout << "       tube radius R = " << R_tube << ", tube separation = " << tube_sep << endl;

    // ======================================================================
    // Build GeneralPatch with TWO types
    // ======================================================================

    vector1<int> vec1(2);
    vec1[0] = 5;        // patches per tube particle
    vec1[1] = M_arms;   // patches per star polymer

    // num_per_type[t] is the (exclusive) upper index of type t
    vector1<int> numb(2);
    numb[0] = N_tube;            // tube particles: 0 .. N_tube-1
    numb[1] = N_tube + N_star;   // star polymers: N_tube .. N_tube+N_star-1

    // -----------------------------
    // Body-frame patch directions
    //   rows 0..4         : tube patches (type 0)
    //   rows 5..5+M_arms-1: star patches (type 1)
    // -----------------------------
    const double lat_x = -sin(pid / N_ring);
    const double lat_y =  cos(pid / N_ring);

    matrix<double> patch_orient(5 + M_arms, 3);

    // Tube patches
    patch_orient(0, 0) = 0.0;    patch_orient(0, 1) = 0.0;     patch_orient(0, 2) =  1.0; // north
    patch_orient(1, 0) = 0.0;    patch_orient(1, 1) = 0.0;     patch_orient(1, 2) = -1.0; // south
    patch_orient(2, 0) = lat_x;  patch_orient(2, 1) =  lat_y;  patch_orient(2, 2) =  0.0; // left
    patch_orient(3, 0) = lat_x;  patch_orient(3, 1) = -lat_y;  patch_orient(3, 2) =  0.0; // right
    patch_orient(4, 0) = 1.0;    patch_orient(4, 1) =  0.0;    patch_orient(4, 2) =  0.0; // binding site (radial outward)

    // Star polymer arm-end patches: Fibonacci spiral on the unit sphere
    if (M_arms == 1) {
        patch_orient(5, 0) = 1.0; patch_orient(5, 1) = 0.0; patch_orient(5, 2) = 0.0;
    } else {
        const double golden = pid * (3.0 - sqrt(5.0));
        for (int a = 0; a < M_arms; a++) {
            double y = 1.0 - (2.0 * a) / (double)(M_arms - 1);
            double rr = sqrt(max(0.0, 1.0 - y * y));
            double th = golden * a;
            patch_orient(5 + a, 0) = cos(th) * rr;
            patch_orient(5 + a, 1) = y;
            patch_orient(5 + a, 2) = sin(th) * rr;
        }
    }

    // -----------------------------
    // Interaction matrix
    //   block layout (rows):
    //     [0 .. 24]                 type0-type0    (5*5)
    //     [25 .. 25+5M-1]           type0-type1    (5*M)
    //     [25+5M .. 25+5M+M^2-1]    type1-type1    (M*M)
    //   within type0-type0:   row = a*5 + b
    //   within type0-type1:   row = 25 + a*M + b   (a in type0, b in type1)
    // -----------------------------
    const int n_pots = 25 + 5 * M_arms + M_arms * M_arms;
    matrix<double> params(n_pots, 3);
    for (int i = 0; i < n_pots; i++) {
        params(i, 0) = 0.0;
        params(i, 1) = range;
        params(i, 2) = patchsize;
    }

    // tube-tube axial & lateral bonds
    params(0 * 5 + 1, 0) = axial_strength;
    params(1 * 5 + 0, 0) = axial_strength;
    params(2 * 5 + 3, 0) = lateral_strength;
    params(3 * 5 + 2, 0) = lateral_strength;

    // tube binding site (patch 4 of type 0) with every arm-end patch of type 1
    for (int b = 0; b < M_arms; b++) {
        params(25 + 4 * M_arms + b, 0) = binding_strength;
    }
    // star-star: leave all zeros (no patchy interaction between stars)

    GeneralPatch c(vec1, numb, params, patch_orient);

    // ======================================================================
    // Initial positions and orientations
    // ======================================================================

    matrix<double> pos(N, 3);
    matrix<double> moms(N, 3);
    matrix<double> orients(N, 9);
    matrix<double> angmoms(N, 3);

    const double z0       = L / 2.0 - (N_z - 1) * dz / 2.0;
    const double y_center = L / 2.0;
    const double x_center = L / 2.0;

    // Two tubes side by side along the x axis
    double tube_x[2];
    tube_x[0] = x_center - tube_sep / 2.0;
    tube_x[1] = x_center + tube_sep / 2.0;

    for (int t = 0; t < N_tubes; t++) {
        for (int j = 0; j < N_z; j++) {
            for (int i = 0; i < N_ring; i++) {
                int k = t * N_tube_each + j * N_ring + i;
                double phi = 2.0 * pid * i / N_ring;

                pos(k, 0) = tube_x[t] + R_tube * cos(phi);
                pos(k, 1) = y_center + R_tube * sin(phi);
                pos(k, 2) = z0 + j * dz;

                // Orientation: body-x aligned with the outward radial direction
                orients(k, 0) =  cos(phi);  orients(k, 1) = sin(phi);  orients(k, 2) = 0.0;
                orients(k, 3) = -sin(phi);  orients(k, 4) = cos(phi);  orients(k, 5) = 0.0;
                orients(k, 6) =  0.0;       orients(k, 7) = 0.0;        orients(k, 8) = 1.0;
            }
        }
    }

    // Star polymers: random positions + random orientations
    // (random-orientation construction copied from LangevinNVTR::initialize)
    for (int s = 0; s < N_star; s++) {
        int k = N_tube + s;

        pos(k, 0) = L * (double)rand() / (double)RAND_MAX;
        pos(k, 1) = L * (double)rand() / (double)RAND_MAX;
        pos(k, 2) = L * (double)rand() / (double)RAND_MAX;

        double x1 = (double)rand() / (double)RAND_MAX;
        double x2 = (double)rand() / (double)RAND_MAX;
        double x3 = (double)rand() / (double)RAND_MAX;
        double v1 = cos(2.0 * pid * x2) * sqrt(x3);
        double v2 = sin(2.0 * pid * x2) * sqrt(x3);
        double v3 = sqrt(1.0 - x3);
        double r1 = cos(2.0 * pid * x1);
        double r2 = sin(2.0 * pid * x1);

        orients(k, 0) = -(r1 * (1 - 2 * SQR(v1))) - 2 * r2 * v1 * v2;
        orients(k, 1) = -(r2 * (1 - 2 * SQR(v1))) + 2 * r1 * v1 * v2;
        orients(k, 2) =  2 * v1 * v3;
        orients(k, 3) =  2 * r1 * v1 * v2 + r2 * (1 - 2 * SQR(v2));
        orients(k, 4) =  2 * r2 * v1 * v2 - r1 * (1 - 2 * SQR(v2));
        orients(k, 5) =  2 * v2 * v3;
        orients(k, 6) =  2 * r1 * v1 * v3 - 2 * r2 * v2 * v3;
        orients(k, 7) =  2 * r2 * v1 * v3 + 2 * r1 * v2 * v3;
        orients(k, 8) = -1 + 2 * SQR(v3);
    }

    // ======================================================================
    // Langevin integrator
    // ======================================================================

    vector1<bool> pb(3, true);
    cube geo(L, pb, 3);

    LangevinNVTR obj(geo);
    obj.initialize(pos, moms, orients, angmoms);
    obj.setdt(dt);
    obj.setgamma(6.0 * pid * viscosity * hdradius);
    obj.setgammar(8.0 * pid * viscosity * hdradius * hdradius * hdradius);
    obj.setkT(kT);

    // ======================================================================
    // Run loop
    // ======================================================================

    int num_cells = (int)floor(L / 4.0);
    int ccc;
    matrix<int> boxes = geo.generate_boxes_relationships(num_cells, ccc);
    matrix<int> *pairs = obj.calculatepairs(boxes, 3.5);

    WCAPotential wsa(3.0, 1.0, 0.0);

    matrix<double> F(N, 3);
    matrix<double> T(N, 3);
    matrix<double> RT(N, 6);

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
