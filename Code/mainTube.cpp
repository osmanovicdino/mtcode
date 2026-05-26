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
//  Two hollow tubes + binding-site lattice + bead-spring star polymers
// ==========================================================================
//
//  TUBE PARTICLES (type 0) — 5 patches each
//    patch 0:  north (+z body axis)       -> axial bond to ring above
//    patch 1:  south (-z body axis)       -> axial bond to ring below
//    patch 2:  lateral "left"             -> bond to CCW neighbour
//    patch 3:  lateral "right"            -> bond to CW  neighbour
//    patch 4:  radial outward (body +x)   -> binding site for star arms
//
//  INNER STAR BEADS (type 1) — 1 patch each (strength = 0, structurally inert)
//    Includes the central bead and all non-terminal arm beads.
//
//  ARM-END STAR BEADS (type 2) — 1 patch each (body +x, binding_strength)
//    Terminal bead of each arm; can bind to tube binding sites (patch 4).
//
//  STAR POLYMER TOPOLOGY
//    Each star: 1 central bead + M_arms arms of length arm_len beads each.
//      total beads per star = 1 + M_arms * arm_len
//      inner beads per star = 1 + M_arms * (arm_len - 1)   (type 1)
//      end   beads per star = M_arms                        (type 2)
//    FENE bonds connect:
//      central -> first bead of each arm, then sequential along each arm.
//
//  INTERACTION TABLE (only non-zero entries):
//    tube  patch 0  <->  tube  patch 1    axial bond
//    tube  patch 2  <->  tube  patch 3    lateral bond
//    tube  patch 4  <->  arm-end patch 0  binding to star
//    all inner star patches: strength 0 (no patchy interaction)
// ==========================================================================

int main(int argc, char **argv)
{
    srand(time(NULL));
    signal(SIGSEGV, handler);

    // ----------------------------------------------------------------------
    // Tube geometry
    // ----------------------------------------------------------------------
    const int    N_ring      = 12;
    const int    N_z         = 60;
    const int    N_tubes     = 2;
    const int    N_tube_each = N_ring * N_z;
    const int    N_tube      = N_tubes * N_tube_each;
    const double chord       = 1.1;
    const double R_tube      = chord / (2.0 * sin(pid / N_ring));
    const double dz          = chord;

    // ----------------------------------------------------------------------
    // Star polymer topology
    // ----------------------------------------------------------------------
    const int    N_star      = 20;
    const int    M_arms      = 6;
    const int    arm_len     = 3;   // beads per arm (includes end bead)
    const int    beads_per_star    = 1 + M_arms * arm_len;
    const int    inner_per_star    = 1 + M_arms * (arm_len - 1);
    const int    end_per_star      = M_arms;
    const int    N_inner     = N_star * inner_per_star;
    const int    N_end       = N_star * end_per_star;
    const int    N           = N_tube + N_inner + N_end;

    // Strides for indexing into the inner-bead and end-bead blocks
    // Central bead of star s:  idx_inner(s, 0)
    // Inner arm bead k of arm m of star s:  idx_inner(s, 1 + m*(arm_len-1) + k)
    // End bead of arm m of star s:  idx_end(s, m)
    auto idx_inner = [&](int s, int local) { return N_tube + s * inner_per_star + local; };
    auto idx_end   = [&](int s, int m)     { return N_tube + N_inner + s * end_per_star + m; };

    // ----------------------------------------------------------------------
    // FENE bond list
    // ----------------------------------------------------------------------
    const int total_fene_bonds = N_star * M_arms * arm_len;
    matrix<int> star_bonds(total_fene_bonds, 2);
    {
        int bi = 0;
        for (int s = 0; s < N_star; s++) {
            int center = idx_inner(s, 0);
            for (int m = 0; m < M_arms; m++) {
                // center -> first inner bead of arm m (or directly to end if arm_len==1)
                int first_in_arm = (arm_len > 1) ? idx_inner(s, 1 + m * (arm_len - 1))
                                                  : idx_end(s, m);
                star_bonds(bi, 0) = center;
                star_bonds(bi, 1) = first_in_arm;
                bi++;

                // sequential bonds along the arm
                for (int k = 0; k < arm_len - 1; k++) {
                    int from = (k < arm_len - 2) ? idx_inner(s, 1 + m * (arm_len - 1) + k)
                                                  : idx_inner(s, 1 + m * (arm_len - 1) + k);
                    // last inner bead -> end bead
                    int to   = (k < arm_len - 2) ? idx_inner(s, 1 + m * (arm_len - 1) + k + 1)
                                                  : idx_end(s, m);
                    star_bonds(bi, 0) = from;
                    star_bonds(bi, 1) = to;
                    bi++;
                }
            }
        }
    }

    // ----------------------------------------------------------------------
    // Box and dynamics
    // ----------------------------------------------------------------------
    const double tube_sep = 3.0 * R_tube + 4.0;
    const double L        = 100.0;
    const double kT       = 1.0;
    const double dt       = 0.005;
    const double viscosity= 1.0;
    const double hdradius = 0.5;

    // ----------------------------------------------------------------------
    // Interaction strengths
    // ----------------------------------------------------------------------
    const double axial_strength   = 15.0;
    const double lateral_strength = 15.0;
    const double binding_strength = 8.0;
    const double range            = 1.2;
    const double patchsize        = 0.6;

    cout << "Setup: " << N_tubes << " tubes x " << N_tube_each
         << " tube particles + " << N_star << " star polymers ("
         << M_arms << " arms x " << arm_len << " beads) = "
         << N << " total particles" << endl;
    cout << "       " << total_fene_bonds << " FENE bonds" << endl;
    cout << "       tube radius R = " << R_tube << endl;

    // ======================================================================
    // GeneralPatch — three types
    //   type 0: tube particles     (5 patches)
    //   type 1: inner star beads   (1 patch, strength 0)
    //   type 2: arm-end star beads (1 patch, binding_strength)
    //
    // Params matrix layout (blocks in upper-triangular order):
    //   (0,0): rows   0..24   tube-tube          5*5 = 25
    //   (0,1): rows  25..29   tube-inner         5*1 =  5
    //   (0,2): rows  30..34   tube-end           5*1 =  5
    //   (1,1): row   35       inner-inner        1*1 =  1
    //   (1,2): row   36       inner-end          1*1 =  1
    //   (2,2): row   37       end-end            1*1 =  1
    //   total: 38 rows
    // ======================================================================

    vector1<int> vec1(3);
    vec1[0] = 5;  vec1[1] = 1;  vec1[2] = 1;

    vector1<int> numb(3);
    numb[0] = N_tube;
    numb[1] = N_tube + N_inner;
    numb[2] = N_tube + N_inner + N_end;

    // Body-frame patch directions
    //   row 0..4: tube patches (north, south, left, right, binding)
    //   row 5:    inner bead patch (dummy direction, strength=0)
    //   row 6:    end bead patch (body +x, faces outward for binding)
    const double lat_x = -sin(pid / N_ring);
    const double lat_y =  cos(pid / N_ring);

    matrix<double> patch_orient(7, 3);
    patch_orient(0, 2) =  1.0;                                        // north
    patch_orient(1, 2) = -1.0;                                        // south
    patch_orient(2, 0) = lat_x;  patch_orient(2, 1) =  lat_y;        // left
    patch_orient(3, 0) = lat_x;  patch_orient(3, 1) = -lat_y;        // right
    patch_orient(4, 0) = 1.0;                                         // binding site
    patch_orient(5, 0) = 1.0;                                         // inner dummy
    patch_orient(6, 0) = 1.0;                                         // end-bead (binding)

    const int n_pots = 38;
    matrix<double> params(n_pots, 3);
    for (int i = 0; i < n_pots; i++) {
        params(i, 0) = 0.0;
        params(i, 1) = range;
        params(i, 2) = patchsize;
    }

    // tube-tube axial bonds: patch 0 (north) <-> patch 1 (south)
    params(0 * 5 + 1, 0) = axial_strength;
    params(1 * 5 + 0, 0) = axial_strength;

    // tube-tube lateral bonds: patch 2 (left) <-> patch 3 (right)
    params(2 * 5 + 3, 0) = lateral_strength;
    params(3 * 5 + 2, 0) = lateral_strength;

    // tube binding site (patch 4 of type 0) <-> arm-end patch (patch 0 of type 2)
    // block (0,2) starts at row 30, patch 4 of type 0 is at offset 4*1+0 = 4
    params(30 + 4, 0) = binding_strength;

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

    double tube_x[2];
    tube_x[0] = x_center - tube_sep / 2.0;
    tube_x[1] = x_center + tube_sep / 2.0;

    // Tube particles: orient body-x radially outward (Rz(-phi))
    for (int t = 0; t < N_tubes; t++) {
        for (int j = 0; j < N_z; j++) {
            for (int i = 0; i < N_ring; i++) {
                int k   = t * N_tube_each + j * N_ring + i;
                double phi = 2.0 * pid * i / N_ring;

                pos(k, 0) = tube_x[t] + R_tube * cos(phi);
                pos(k, 1) = y_center  + R_tube * sin(phi);
                pos(k, 2) = z0 + j * dz;

                orients(k, 0) =  cos(phi);  orients(k, 1) = sin(phi);  orients(k, 2) = 0.0;
                orients(k, 3) = -sin(phi);  orients(k, 4) = cos(phi);  orients(k, 5) = 0.0;
                orients(k, 6) =  0.0;       orients(k, 7) = 0.0;       orients(k, 8) = 1.0;
            }
        }
    }

    // Star polymers: central bead placed randomly, arms along Fibonacci-spiral directions,
    // beads spaced by bond_length so that FENE bonds start near equilibrium.
    const double bond_length = 0.9;
    const double golden_star = pid * (3.0 - sqrt(5.0));

    for (int s = 0; s < N_star; s++) {
        double cx = L * (double)rand() / RAND_MAX;
        double cy = L * (double)rand() / RAND_MAX;
        double cz = L * (double)rand() / RAND_MAX;

        // Central bead (type 1, inner) — identity orientation
        int ci = idx_inner(s, 0);
        pos(ci, 0) = cx;  pos(ci, 1) = cy;  pos(ci, 2) = cz;
        orients(ci, 0) = 1.0; orients(ci, 4) = 1.0; orients(ci, 8) = 1.0;

        for (int m = 0; m < M_arms; m++) {
            // Arm direction from Fibonacci spiral
            double yy = (M_arms > 1) ? 1.0 - (2.0 * m) / (double)(M_arms - 1) : 0.0;
            double rr = sqrt(max(0.0, 1.0 - yy * yy));
            double th = golden_star * m;
            double adx = cos(th) * rr;
            double ady = yy;
            double adz = sin(th) * rr;

            // Inner arm beads (arm_len - 1 beads)
            for (int k = 0; k < arm_len - 1; k++) {
                int ii = idx_inner(s, 1 + m * (arm_len - 1) + k);
                pos(ii, 0) = cx + (k + 1) * bond_length * adx;
                pos(ii, 1) = cy + (k + 1) * bond_length * ady;
                pos(ii, 2) = cz + (k + 1) * bond_length * adz;
                orients(ii, 0) = 1.0; orients(ii, 4) = 1.0; orients(ii, 8) = 1.0;
            }

            // End bead (type 2) — orient body-x along arm direction (outward)
            int ei = idx_end(s, m);
            pos(ei, 0) = cx + arm_len * bond_length * adx;
            pos(ei, 1) = cy + arm_len * bond_length * ady;
            pos(ei, 2) = cz + arm_len * bond_length * adz;

            // Build rotation matrix with body-x = arm direction
            // Construct an orthonormal frame: e1=arm_dir, e2, e3
            double e1x = adx, e1y = ady, e1z = adz;
            // Pick an arbitrary perpendicular vector
            double e2x, e2y, e2z;
            if (fabs(e1x) < 0.9) {
                e2x = 0.0; e2y = -e1z; e2z = e1y;
            } else {
                e2x = -e1z; e2y = 0.0; e2z = e1x;
            }
            double e2n = sqrt(e2x*e2x + e2y*e2y + e2z*e2z);
            e2x /= e2n; e2y /= e2n; e2z /= e2n;
            double e3x = e1y*e2z - e1z*e2y;
            double e3y = e1z*e2x - e1x*e2z;
            double e3z = e1x*e2y - e1y*e2x;

            // R stored row-major: row0=e1, row1=e2, row2=e3
            orients(ei, 0) = e1x; orients(ei, 1) = e1y; orients(ei, 2) = e1z;
            orients(ei, 3) = e2x; orients(ei, 4) = e2y; orients(ei, 5) = e2z;
            orients(ei, 6) = e3x; orients(ei, 7) = e3y; orients(ei, 8) = e3z;
        }
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
    // Potentials
    // ======================================================================

    WCAPotential  wca(3.0, 1.0, 0.0);
    FENEPotential fene(30.0, 1.5);  // kbond=30, R_0=1.5

    // ======================================================================
    // Run loop
    // ======================================================================

    outfunc(obj.getdat(),"pos.csv");
    pausel();

    const int runtime = 10000000;
    const int every   = 1000;

    int num_cells = (int)floor(L / 4.0);
    int ccc;
    matrix<int>  boxes = geo.generate_boxes_relationships(num_cells, ccc);
    matrix<int> *pairs = obj.calculatepairs(boxes, 3.5);

    matrix<double> F(N, 3);
    matrix<double> T(N, 3);
    matrix<double> RT(N, 6);

    // Initial force/torque evaluation
    F  = obj.calculateforces(*pairs, wca);
    F += obj.calculateforces(star_bonds, fene);
    obj.calculate_forces_and_torques3D(*pairs, c, F, T);
    generate_uniform_random_matrix(RT);
    obj.create_forces_and_torques_sphere(F, T, RT);

    int number_of_digits = 0;
    {
        int tf = (runtime / every) + 1;
        do { ++number_of_digits; tf /= 10; } while (tf);
    }

    for (int i = 0; i < runtime; i++) {
        if (i % every == 0)
        {
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << (i / every);
            outfunc(obj.getdat(), "tube_pos_i=" + ss.str() );
            outfunc(obj.getorientation(), "tube_ori_i=" + ss.str() );
        }

        if (i > 0 && i % 20 == 0) {
            delete pairs;
            pairs = obj.calculatepairs(boxes, 3.5);
        }

        obj.advancemom_halfstep(F, T);
        obj.advance_pos();
        obj.rotate();

        F  = obj.calculateforces(*pairs, wca);
        F += obj.calculateforces(star_bonds, fene);
        T.reset(0.0);
        obj.calculate_forces_and_torques3D(*pairs, c, F, T);
        generate_uniform_random_matrix(RT);
        obj.create_forces_and_torques_sphere(F, T, RT);
        obj.advancemom_halfstep(F, T);


    }

    delete pairs;
    return 0;
}