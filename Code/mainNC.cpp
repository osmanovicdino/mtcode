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

#include "basic.h"
#include "vector1.h"
#include "matrix2.h"
#include "matrix2.cpp"
#include "potential.h"
//#include "intmatrix.h"
#include "MD.h"
#include "Langevin.h"


// #include "BrownianGel.cpp"
// #include "BrownianGel2.cpp"
// #include "LangevinGel.cpp"
// #include "LangevinGelFixed.cpp"

#include "NCGasR.h"



using namespace std;




vector1<int> calcg1(matrix<double> &data, cube &geo, double max, double dx) { 

vector<double> index1;


double ll  = geo.l;
int ccc;
int num = floor(ll/max); //NB THIS ALGORITHM GOES WRONG IF THERE ARE LESS THAN 3 BOXES num=2

matrix<int> boxlist = geo.generate_boxes_relationships(num,ccc);

int total_cubes = boxlist.getNsafe();

double dims = (double)(geo.dimension);

int cubes_per_length = (int)round(exp(log(total_cubes)/dims));


vector<vector<int> > b;
b.reserve(total_cubes);
for(int j = 0 ; j < total_cubes ; j++) {
    vector<int> temp;
    b.push_back(temp);
}
int dimension = (geo.dimension);


vector1<int> dim(dimension);
for(int i = 0 ; i < dimension ; i++ ) {
	int ij = 1;
	for(int j = 0 ; j < i ; j++ ) {
	ij*= cubes_per_length;
	}
	dim[i] = ij;
}




for(int i = 0 ; i < data.getNsafe() ; i++) {



	int c = geo.assign_box(data,i,dim,cubes_per_length);

	b[c].push_back(i);
}



int ss = boxlist.getncols();



	#pragma omp parallel
	{
	vector<double> index1_private;

	#pragma omp for nowait schedule(static)
	for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++)
		for(int c2 = 0 ; c2 < ss ; c2++) {

			int box1 = c1;
			int box2 = boxlist(c1,c2);
			if(box1==box2) {
				for(int i = 0 ; i < b[box1].size() ; i++) {
					for(int j = i+1 ;  j < b[box2].size() ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];
						//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);

							//if(int_dl[q]) {
						double dis =  geo.distance(data,iterator1,iterator2);
		
						
							
						//int q = ints.get_potential_number(iterator1,iterator2,l);
						index1_private.push_back(dis);
						//index2_private.push_back(iterator2);
						//index3_private.push_back(l);	
											
													
					}
				}
			}
			else if(box2>box1) {
				for(int i = 0 ; i < b[box1].size() ; i++) {
					for(int j = 0 ;  j < b[box2].size() ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];

						double dis =  geo.distance(data,iterator1,iterator2);
							
						//int q = ints.get_potential_number(iterator1,iterator2,l);
						index1_private.push_back(dis);					
						// for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
						// 	int q = ints.get_potential_number(iterator1,iterator2,l);

						// 	if(int_dl[q]) {
						// 		bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
						// 		if(cond) {
						// 			index1_private.push_back(iterator1);
						// 			index2_private.push_back(iterator2);
						// 			index3_private.push_back(l);					
						// 		}

						// 	}
						// 	else {
						// 		index1_private.push_back(iterator1);
						// 		index2_private.push_back(iterator2);
						// 		index3_private.push_back(l);
						// 	}							
						// }						
					}
				}
			}
			else {

			}
	}


	#pragma omp for schedule(static) ordered
	for(int i = 0 ; i < omp_get_num_threads(); i++) {
		#pragma omp ordered
		index1.insert(index1.end(),index1_private.begin(),index1_private.end());
	}

	}

	vector1<double> a(index1.size());
	for(int i = 0 ; i < a.getsize() ; i++) {
		a[i] = index1[i];
	}
	//cout << "new" << endl;
	//create  a histogram of all the distances
	int gs = floor(max/dx);
	vector1<int> g(gs);
	//where index i is (i)*dx to (i+1)*dx

	for(int i = 0  ; i < a.getsize() ; i++) {
		int in = floor(a[i]/dx);
		if(in<gs) g[in]+=1;
	}	


	return g;


}


vector1<double> densitycorrelations(matrix<double> &data, cube &geo, double max, int howm) {
vector<double> index1;

double ll  = geo.l;
int ccc;
// int num = floor(ll/mg);
int num  = howm;

double dx = ll/(double)howm;

matrix<int> boxlist = geo.generate_boxes_relationships(num,ccc);

int total_cubes = boxlist.getNsafe();

double dims = (double)(geo.dimension);

int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

vector<vector<int> > b;
b.reserve(total_cubes);
for(int j = 0 ; j < total_cubes ; j++) {
    vector<int> temp;
    b.push_back(temp);
}
int dimension = (geo.dimension);


vector1<int> dim(dimension);
for(int i = 0 ; i < dimension ; i++ ) {
	int ij = 1;
	for(int j = 0 ; j < i ; j++ ) {
	ij*= cubes_per_length;
	}
	dim[i] = ij;
}



for(int i = 0 ; i < data.getNsafe() ; i++) {
	int c = geo.assign_box(data,i,dim,cubes_per_length);
	b[c].push_back(i);
}

matrix<double> as(b.size(),4);


for(int i = 0 ; i < cubes_per_length; i++ )
	for(int j = 0; j < cubes_per_length; j++)
		for(int k = 0 ; k < cubes_per_length ; k++) {
			vector1<int> scal(3);
			int numc = i+cubes_per_length*j+cubes_per_length*cubes_per_length*k;
			as(numc,0)=0.5*dx+i*dx;
			as(numc,1)=0.5*dx+j*dx;
			as(numc,2)=0.5*dx+k*dx;
			as(numc,3)=b[numc].size();
		}


	int tot =  30;

	double lk = 0.86602540378*ll;

	double dx2 = lk/(double)tot;



	vector1<double> rhorho(tot);
	vector1<double> val(tot);

	for(int i = 0 ; i < as.getNsafe() ; i++ ) {
		for(int j = i ; j < as.getNsafe() ; j++) {
			double dis =  geo.distance(as,i,j);
			int in = floor(dis/dx2);
			if(in<tot) {
				val[in]+=1;
				rhorho[in]+=as(i,3)*as(j,3);
			}	
		}
	}

	return rhorho/val;

// cout << b.size() << endl;
// cout << b[0].size() << endl;
// pausel();



}



int main(int argc, char** argv) {
srand (time(NULL));

// double eps[3] = {0.0,1.0,2.0};
// double eqeps[3] = {0.0,1.0,2.0};
// double dis[10] = {29.6931, 21.8781, 17.3647, 15.1694, 13.7823, 12.7944, 12.04, 11.4369, 10.939, 10.1549};
		// NCGas a(80.572,10000);

		// double c1 = 0.1;
		// double c2 = 2.0;

		// cout << c1 << ", " << c2 << endl;

		// a.seteps(c1);

		// a.setkT(1.0);

		// a.seteqeps(c2);

		// a.run(2000000);

	int n = 10000;
	double v00 =  0.0;
	//for(double v00 = 0.0 ; v00 < 101.0 ; v00 += 25.0) {
	for(double eqeps =  2.1 ; eqeps < 2.35 ; eqeps += 0.1 ) {
	NCGasR a(140.0,n,2);

	a.seteqeps(eqeps);
	a.setv00(v00);

	cout << "here" << endl;

	a.run2D(10000000);
	}
//	}
	
// ofstream myfile;
// myfile.open("/home/dino/Code/DinoFastMD/data_from_sims2/paircorrelations.txt");
			
// // double eps[3] = {0.0,1.5,3.};
// // double eqeps[3] = {0.0,0.5,1.0};
// // //double kT[4] = {0.5,1.0,1.5,2.0};
// vector1<bool> pb(3,true);
// // // NCGas a(12.04,1000);
// NCGas a(10.0,500);

// for(int i = 0 ; i < 3 ; i++) {
// 	for(int j = 0  ; j < 3 ; j++ ) {
		

// 		for(int k = 0 ; k < 10 ; k++ ) {
// 			double leng = dis[k];
			
// 			double c1 = eps[i];
// 			double c2 = eqeps[j];
// 			double c3 = 1.0;

// 			//cout << c1 << ", " << c2 << endl;
// 			vector1<bool> pb(3,true);
// 			cube bc(leng,pb,3);

// 			a.setgeo(bc);

// 			a.seteps(c1);

// 			a.seteqeps(c2);

// 			a.setkT(c3);

// 			for(int it = 1 ; it < 2000 ; it++) {
			
			


// 			stringstream kts;
// 		 	kts << (*(a.obj)).getkT();
			
// 			 stringstream epi;
// 			 epi << c1;

// 			 stringstream epieq;
// 			 epieq << c2;

// 			 stringstream len;
// 			 len << leng;

// 			 string extension =  "_eqeps="+epieq.str()+"_eps="+epi.str()+"_kT="+kts.str()+"_l="+len.str()+".csv";
// 			// pairlist += ss2.str();
// 			// pairlist += extension;

// 			// ofstream myfilex;
// 			// myfilex.open(pairlist.c_str());
// 			// myfilex <<= froyo;
// 			// myfilex.close();		
		
// 			stringstream ss;
// 			ss <<it;
// 			string filename = "x";
// 			filename += ss.str();
// 			filename += extension;

// 			string forcename = "F";
// 			forcename += ss.str();
// 			forcename += extension;

// 			string fotalfilename = "/home/dino/Code/DinoFastMD/data_from_sims2/"+forcename;

// 			string totalfilename = "/home/dino/Code/DinoFastMD/data_from_sims2/"+filename;

// 			cout << totalfilename << endl;
// 			cout << fotalfilename << endl;


// 			double T;
// 			bool err1;
// 			matrix<double> state = importcsv(totalfilename,T,err1);
// 			bool err2;
// 			matrix<double> ff =  importcsv(fotalfilename,T,err2);


// 			if(err1) {cout << "read error" << endl; exit(1); }
// 			if(err2) {cout << "read error" << endl; exit(1); }

// 			(*(a.obj)).setdat(state);
// 			cout << "state set" << endl;
// 			//cube bc(leng,pb,3);
// 			vector1<int> g1 = calcg1(state, bc, 3.0, 0.05);
// 			cout << "g1" << endl;
// 			matrix<double> sigma(3,3);// = a.calculate_virial();


// 			myfile << totalfilename << "\t";
// 			myfile << sigma(0,0) << "," << sigma(0,1) << "," << sigma(0,2) << "," << sigma(1,0) << "," << sigma(1,1) << "," << sigma(1,2) << "," << sigma(2,0) << "," << sigma(2,1) << "," << sigma(2,2) << "\t";
// 			myfile <<= g1;
// 			myfile << endl;
			



// 			}
// 		}
// 	}
// }
// myfile.close();


		//double len = 10.939;
		
		// NCGas a(len,500);

		// double c1 = 2.0;
		// double c2 = 2.0;

		// cout << c1 << ", " << c2 << endl;

		// a.seteps(c1);

		// a.setkT(1.0);

		// a.seteqeps(c2);

		// a.run(100000);


// double eps[3] = {0.0,1.0,2.0};
// double eqeps[3] = {0.0,1.0,2.0};
// double dis[10] = {29.6931, 21.8781, 17.3647, 15.1694, 13.7823, 12.7944, 12.04, 11.4369, 10.939, 10.1549};
// vector<vector1<double> > datapoints;
// datapoints.reserve(10*3*3);
// for(int i = 0  ; i < 3 ; i++) {
// 	for(int j = 0  ; j < 3 ; j++) {
// 		for(int k = 0 ; k < 10 ; k++ ) {
// 			vector1<double> temp(3);
// 			temp[0]=eps[i];
// 			temp[1]=eqeps[j];
// 			temp[2]=dis[k];
// 			datapoints.push_back(temp);
// 		}
// 	}
// }



// int simrun = 5;
// int start_i = simrun*15;
// int end_i = start_i+15;

// for(int i = start_i ; i < end_i ; i++ ) {
// 		vector1<double> param = datapoints[i];

// 		double len = param[2];
// 		NCGas a(len,500);

// 		double c1 = param[0];
// 		double c2 = param[1];

// 		cout << c1 << ", " << c2 << endl;

// 		a.seteps(c1);

// 		a.setkT(1.0);

// 		a.seteqeps(c2);

// 		a.run(2000000);
// 	}

// int start_i=1;
// int end_i=start_i+1;
// int start_k=0;
// int end_k=5;



// for(int i = start_i ; i < end_i ; i++ ) {
// 	for(int j = 0 ; j < 1 ; j++) {
// 		for(int k = start_k  ; k < end_k; k++) {
// 		double len = dis[k];
// 		NCGas a(len,500);		

// 		double c1 = eps[i];
// 		double c2 = eqeps[j];

// 		cout << c1 << ", " << c2 << endl;

// 		a.seteps(c1);

// 		a.setkT(1.0);

// 		a.seteqeps(c2);

// 		a.run(3000000);
// 		}
// 	}
// }


// double eps[4] = {0.0,1.0,2.0,3.0};
// double eqeps[1] = {0.0};
// double dis[4] = {12.04, 11.4369, 10.939, 10.1549};
// vector<vector1<double> > datapoints;
// datapoints.reserve(4*1*4);
// for(int i = 0  ; i < 4 ; i++) {
// 	for(int j = 0  ; j < 1 ; j++) {
// 		for(int k = 0 ; k < 4 ; k++ ) {
// 			vector1<double> temp(3);
// 			temp[0]=eps[i];
// 			temp[1]=eqeps[j];
// 			temp[2]=dis[k];
// 			datapoints.push_back(temp);
// 		}
// 	}
// }



// int simrun = 3;
// int start_i = simrun*4;
// int end_i = start_i+4;

// for(int i = start_i ; i < end_i ; i++ ) {
// 		vector1<double> param = datapoints[i];

// 		double len = param[2];
// 		NCGas a(len,500);

// 		double c1 = param[0];
// 		double c2 = param[1];

// 		cout << c1 << ", " << c2 << endl;

// 		a.seteps(c1);

// 		a.setkT(1.0);

// 		a.seteqeps(c2);

// 		a.run(2000000);
// 	}

return 0;
}
