#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/ioctl.h> 
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <random>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "Basic/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"
#include "MDBase/Potentials/potential.h"
//#include "intmatrix.h"
#include "MDBase/MD.h"
#include "MDBase/Langevin.h"


// #include "BrownianGel.cpp"
// #include "BrownianGel2.cpp"
// #include "LangevinGel.cpp"
// #include "LangevinGelFixed.cpp"

// #include "NCGasR.h"
#include "Microtubule/Microtubule.h"


using namespace std;






int main(int argc, char** argv) {
srand (time(NULL));



//SET LENGTH OF THE SIMULATION BOX

double simulation_box_length =  200.0;

//SET THE TOTAL NUMBER OF PARTICLES OF EACH TYPE

int particles_of_typeA = 0;

int particles_of_typeB = 0; //set to 0 for 1 liquid

//SET NUMBER OF MICROTUBULES

string filename1 = "/home/dino/Desktop/tylercollab/LegacyCode/config_files/pa.csv";  //FULL FILE PATH TO INITIAL A PARTICLE DATA
string filename2 = "/home/dino/Desktop/tylercollab/LegacyCode/config_files/pb.csv"; //FULL FILE PATH TO INITIAL B PARTICLE DATA
string filename3 = "/home/dino/Desktop/tylercollab/LegacyCode/config_files/pc.csv"; //FULL FILE PATH TO INITIAL MIRCROTUBULE DATA

int number_of_mircrotubules = 1;

//SET LENGTH OF MICROTUBULE IN MONOMERS

// int length_of_microtubule = 150; //SET THIS TO BE THE TOTAL NUMBER OF PARTICLES IN THE SYSTEM
bool err1;
double tg;
int length_of_microtubule = (importcsv(filename3, tg, err1)).getnrows();
//Call the microtubule constructor;

Microtubule a(simulation_box_length,particles_of_typeA,particles_of_typeB,number_of_mircrotubules,length_of_microtubule);


//in order to change any of the potentials between subparts of the system we can set a new WCA potential between the different components


//WCAPotential potential1(epsilon1,sigma,epsilon2); //this is a WCA potential where epsilon2 is the attractive part of the interaction


// we can change potentials from the following list, uncomment and add whichever potential you want in order to change interactions between parts of the system

// a.setpotaa(potential1); // interaction between A and A (default is WCA with att epsilon =2)
// a.setpotab(potential1); //interaction between A and B (red/blue) (default is steric repulsion)
// a.setpotac(potential1); //interaction between A and microtubule (default is steric repulsion)
// a.setpotbb(potential1); //interaction between B and B (blue/blue) (default is WCA with att epsilon =2)
// a.setpotbc(potential1); //interaction between B and microtubule (default is steric repulsion)
// a.setpotcc(potential1); //interaction between microtubule monomers //default is steric repulsion

//OTHER SIMULATION PARAMETERS

double v00_a = 6.0; // force with which the microtubule is propelled by A particles.
double v00_b = 6.0; // force with which the microtubule is propelled by B particles

a.setv0(v00_a,v00_b);


double pola=1.;
double polb=-1.;

a.setpolarity(pola,polb); //set the polarities of the motors (+1 -> plus end, -1 -> minus end)



//probabilities of binding/unbinding

double probunbind =0.0; //unbinding probability, note that the motor detaches naturally itself when it reaches the ends, so this is an unbinding rate when the motor is not at the ends.
						//if this number is zero this means that the 

double probbind= 1.0; //this sets the rate at which motors within the right distance bind to the microtubule, a rate of one means that whenever a motor is within the correct distance of a microtubule it will bind


double excess_force_distance = 2.; //if the motor and the nanostar are more than this distance apart, the motor will drop off the microtubule

a.setprobunbind(probunbind);
a.setprobbind(probbind);
a.set_excess_force_distance(excess_force_distance);

//in addition to these microtubule parameters, we have other general simulation parameters, such as the viscous damping term:

//UNCOMMENT THE FOLLOWING FOR SETTING UP AN INITIAL SYSTEM


a.set_initial_conditions(filename1,filename2,filename3); //THE ROUTINE TO IMPORT STRING FILENAMES IS IN MATRIX2.CPP (importcsv) if there are architectural differences you need to modify the function

//THE INITIAL DATA SHOULD BE SET UP AS CSV FILE:
/*
particle1x,particle1y
particle2x,particle2y
particle3x,particle3y
...
*/


double gamma = 50.0; //(overdamped)
a.setgamma(gamma);



// struct spatial_viscosity {
// 	double operator()(double,double) { //define the spatially varying intrinsic viscosity
// 		return 50.; //define whichever function you like for the spatial viscosity
// 	}
// };


// spatial_viscosity q;

// bool spatial_variance = false; //false if not spatially varying, true if it is.
//IF FALSE the run will use gamma
//IF TRUE the function will use q (whichever function you define)


// a.set_gamma_spatial_dependence(spatial_variance);


cout << "creation" << endl;

//this will run the simulation, outputting a csv of positions every 1000 timesteps in the directory in which you run the file
 /* UNCOMMENT WHICHEVER SYSTEM YOU WANT, CONSTANT GAMMA, SPATIAL VISOCITY, OR PARTICLE DEPENDENT VISCOSITY */


//CONSTANT GAMMA

int timesteps_between_writing_to_file = 1000;

//a.run(10000000,timesteps_between_writing_to_file);

//a.run(10000000,timesteps_between_writing_to_file);


//SPATIAL VISOCSITY 

//DEFINED AS A INFINITE SPATIAL GRADIENT IN VISCOSITY, above 100.0 it is one value, and less than 100.0 it is another

struct spatial_viscosity {
	double operator()(double x, double y) { //define the spatially varying intrinsic viscosity
		if(x > 100.0 )	return 50.; //define whichever function you like for the spatial viscosity
		else return 25.;
	}
};


// spatial_viscosity q;

// int timesteps_between_writing_to_file = 1000;

// a.runSV(10000000,timesteps_between_writing_to_file,q);

//PARTICLE DEPENDENT VISCOSITY


// struct particle_viscosity {
// 	double operator()(int i) { //define the spatially varying intrinsic viscosity
// 		if(i < 5000 ) {
// 			return 50.0;
// 		}
// 		else {
// 			return 100.0;
// 		}
// 	}
// };


// particle_viscosity q;

// int timesteps_between_writing_to_file = 1000;

// a.runPV(10000000,timesteps_between_writing_to_file,q);

spatial_viscosity q;




cout << "start" << endl;

//This is the distance we create bonds up tp
double distance_to_create_bonds_up_to = 2.0;

//the following vectors define the orientation of the active force:

vector1<int> p1(3);

vector1<int> p2(3);

//these define the directions of the active forces, going along the axis of p1 to p2.
p1[0] = 0; p2[0] = 1;
p1[1] = 2; p2[1] = 3;
p1[2] = 4; p2[2] = 5;


a.runMTONLY_initialstate(10000,1000,q,distance_to_create_bonds_up_to,p1,p2);


return 0;
}
