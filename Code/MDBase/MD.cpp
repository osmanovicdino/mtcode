MD::MD() : geo(cube())
{
	//dat = new matrix<double>();
	//cout << "abstract base class MD called" << endl;
	HSPotential temppot(1.0,0.0);
	ints = temppot.clone();
	dat = new matrix<double>(1,3);
	

}

MD::MD(const MD &old) : geo(old.geo)
{
	// geometry *geotemp = (old.geo);
	// //geo = &(old.geo);
	//cout << "copy constructor MD called" << endl;
	
	// cube *geotemp = (old.geo)->clone();
	// geo = geotemp;

	//cout << "yo" << endl;
	//geo = old.geo;
	//dat = new matrix<double>;
	potential *potnew =  (old.ints)->clone();
	ints = potnew;

	matrix<double> *groleo = (old.dat)->clone();
	dat = groleo;

	dimension = dat->getncols();

	if (dimension != geo.dimension) error("copy constructor dimensions must match in MD");
}

MD& MD::operator=(const MD &old) {
	//cout << "operator= MD called" << endl;

	geo =  old.geo;

	potential *potnew = (old.ints)->clone();
	ints = potnew;

	matrix<double> *groleo = (old.dat)->clone();
	dat = groleo;

	dimension = dat->getncols();

	return *this;
}


MD::~MD() {
	delete dat;
	// delete geo;
}

void MD::setgeometry(cube &a) {
	//cube* q = a.clone();
	geo = a;


}

void MD::setdat(const matrix<double> &a) {
	// matrix<double> *res =  a.clone();
	// dat = res;

	delete dat;

	dat = a.clone();
	dimension = dat->getncols();
	if (dimension != (geo.dimension) ) {
		cout << dimension << " " << geo.dimension << endl;
		error("set dat dimensions must match in MD");
	}

}

void MD::set_particle(const vector1<double> &a, int index)
{
	// matrix<double> *res =  a.clone();
	// dat = res;
	if(a.size() != dimension) error("incorrect dimension in set_particle");
	
	for(int i1 = 0 ; i1 < dimension ; i1++)
		(*dat)(index,i1) = a.gpcons(i1);

}

vector1<double> MD::get_particle(int index)
{
	// matrix<double> *res =  a.clone();
	// dat = res;
	vector1<double> a(dimension);


	for (int i1 = 0; i1 < dimension; i1++)
		a[i1] = (*dat)(index, i1);
	return a;
}

void MD::setinteractions(potential &a) {
	potential* q = a.clone();
	ints = q;
}

double MD::getcoordinate(int i, int j) {
	return (*dat).mat[i*dimension+j];
}

int MD::getdimension() {
	return dimension;
}

int MD::getN() {
	return (*this->dat).getNsafe();
}

matrix<double> MD::getdat() const {
	return (*(this->dat));
}

cube MD::getgeo() const {
	return geo;
}

potential& MD::getints() {
	return *ints;
}

void MD::disvec(int &i1, int &i2, vector1<double> &un, double &dis) {
	geo.distance_vector(*dat,i1,i2,un,dis);
}

double MD::distance(const int &i,const int &j) {
return geo.distance(*dat,i,j);
}

// bool MD::distance_less_than(const int &i,const int &j, double R) {
// return (*this->geo).distance_less_than((*dat)[i],(*dat)[j],R);
// }




// matrix<int>* MD::calculatepairs(matrix<int> &boxlist) {
// 	vector<int> index1;
// 	vector<int> index2;
// 	vector<int> index3;
// 	//#pragma omp parallel for
// 	int Ns = (*dat).getNsafe();
// 	int np = ints.number_of_potentials();
// 	//vector1<double> int_dis(np);
// 	vector1<double> int_dis(np);
// 	vector1<bool> int_dl(np);
	
// 	for(int i = 0 ; i < np ; i++) {
// 		int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
// 	}
	
// 	for(int i = 0 ; i < np ; i++) {
// 		int_dl[i] = (ints.access_potential(i).dl);
// 	}


// 	int total_cubes = boxlist.getNsafe();

// 	double dims = (double)this->getdimension();

// 	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

// 	vector<vector<int> > b;
// 	b.reserve(total_cubes);
//     for(int j = 0 ; j < total_cubes ; j++) {
//         vector<int> temp;
//         b.push_back(temp);
//     }
//     int dimension = this->getdimension();


// 	vector1<int> dim(dimension);
// 	for(int i = 0 ; i < dimension ; i++ ) {
// 		int ij = 1;
// 		for(int j = 0 ; j < i ; j++ ) {
// 		ij*= cubes_per_length;
// 		}
// 		dim[i] = ij;
// 	}



// 	for(int i = 0 ; i < this->getN() ; i++) {

// 		int c = geo.assign_box((*dat)[i],dim,cubes_per_length);

// 		b[c].push_back(i);
// 	}



// 	int ss = boxlist[0].getsize();

// 	#pragma omp parallel
// 	{
// 	vector<int> index1_private;
// 	vector<int> index2_private;
// 	vector<int> index3_private;
// 	#pragma omp for nowait schedule(static)
// 	for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++)
// 		for(int c2 = 0 ; c2 < ss ; c2++) {
// 			int box1 = c1;
// 			int box2 = boxlist(c1,c2);
// 			if(box1==box2) {
// 				for(int i = 0 ; i < b[box1].size() ; i++) {
// 					for(int j = i+1 ;  j < b[box2].size() ; j++) {
// 						int iterator1 = (b[box1])[i];
// 						int iterator2 = (b[box2])[j];
// 						for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							int q = ints.get_potential_number(iterator1,iterator2,l);

// 							if(int_dl[q]) {
// 								bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
// 								if(cond) {
// 									index1_private.push_back(iterator1);
// 									index2_private.push_back(iterator2);
// 									index3_private.push_back(l);					
// 								}

// 							}
// 							else {
// 								index1_private.push_back(iterator1);
// 								index2_private.push_back(iterator2);
// 								index3_private.push_back(l);
// 							}							
// 						}
// 					}
// 				}

// 			}
// 			else if(box2>box1) {
// 				for(int i = 0 ; i < b[box1].size() ; i++) {
// 					for(int j = 0 ;  j < b[box2].size() ; j++) {
// 						int iterator1 = (b[box1])[i];
// 						int iterator2 = (b[box2])[j];
// 						for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							int q = ints.get_potential_number(iterator1,iterator2,l);

// 							if(int_dl[q]) {
// 								bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
// 								if(cond) {
// 									index1_private.push_back(iterator1);
// 									index2_private.push_back(iterator2);
// 									index3_private.push_back(l);					
// 								}

// 							}
// 							else {
// 								index1_private.push_back(iterator1);
// 								index2_private.push_back(iterator2);
// 								index3_private.push_back(l);
// 							}							
// 						}						
// 					}
// 				}
// 			}
// 			else {

// 			}
// 	}
// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index1.insert(index1.end(),index1_private.begin(),index1_private.end());
// 	}
// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index2.insert(index2.end(),index2_private.begin(),index2_private.end());
// 	}

// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index3.insert(index3.end(),index3_private.begin(),index3_private.end());
// 	}
// 	}


// 	matrix<int> * a = new s_matrix<int>(index1.size(),3);
// 	//s_matrix<int> pairs(index1.size(),3);
// 	for(int i = 0 ; i < (*a).getNsafe() ; i++) {
// 		(*a)(i,0) = index1[i];
// 		(*a)(i,1) = index2[i];
// 		(*a)(i,2) = index3[i];
// 	}

	
	
// 	return a;


// }


// matrix<int>* MD::calculatepairs_global(matrix<int> &boxlist, double cut_off) {
// 	vector<int> index1;
// 	vector<int> index2;
// 	vector<int> index3;
// 	//#pragma omp parallel for
// 	// int Ns = (*dat).getNsafe();
// 	// int np = ints.number_of_potentials();
// 	// //vector1<double> int_dis(np);
// 	// vector1<double> int_dis(np);
// 	// vector1<bool> int_dl(np);
	
// 	// for(int i = 0 ; i < np ; i++) {
// 	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
// 	// }
	
// 	// for(int i = 0 ; i < np ; i++) {
// 	// 	int_dl[i] = (ints.access_potential(i).dl);
// 	// }


// 	int total_cubes = boxlist.getNsafe();

// 	double dims = (double)this->getdimension();

// 	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

// 	vector<vector<int> > b;
// 	b.reserve(total_cubes);
//     for(int j = 0 ; j < total_cubes ; j++) {
//         vector<int> temp;
//         b.push_back(temp);
//     }
//     int dimension = this->getdimension();


// 	vector1<int> dim(dimension);
// 	for(int i = 0 ; i < dimension ; i++ ) {
// 		int ij = 1;
// 		for(int j = 0 ; j < i ; j++ ) {
// 		ij*= cubes_per_length;
// 		}
// 		dim[i] = ij;
// 	}



// 	for(int i = 0 ; i < this->getN() ; i++) {

// 		int c = geo.assign_box((*dat)[i],dim,cubes_per_length);

// 		b[c].push_back(i);
// 	}



// 	int ss = boxlist[0].getsize();

// 	#pragma omp parallel
// 	{
// 	vector<int> index1_private;
// 	vector<int> index2_private;
// 	vector<int> index3_private;
// 	#pragma omp for nowait schedule(static)
// 	for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++)
// 		for(int c2 = 0 ; c2 < ss ; c2++) {
// 			int box1 = c1;
// 			int box2 = boxlist(c1,c2);
// 			if(box1==box2) {
// 				for(int i = 0 ; i < b[box1].size() ; i++) {
// 					for(int j = i+1 ;  j < b[box2].size() ; j++) {
// 						int iterator1 = (b[box1])[i];
// 						int iterator2 = (b[box2])[j];
// 						//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							//int q = ints.get_potential_number(iterator1,iterator2,l);

// 							//if(int_dl[q]) {
// 						bool cond =  distance_less_than(iterator1,iterator2,cut_off);
		
// 						if(cond) {
// 							for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							//int q = ints.get_potential_number(iterator1,iterator2,l);
// 							index1_private.push_back(iterator1);
// 							index2_private.push_back(iterator2);
// 							index3_private.push_back(l);	
// 							}				
// 						}							
// 					}
// 				}
// 			}
// 			else if(box2>box1) {
// 				for(int i = 0 ; i < b[box1].size() ; i++) {
// 					for(int j = 0 ;  j < b[box2].size() ; j++) {
// 						int iterator1 = (b[box1])[i];
// 						int iterator2 = (b[box2])[j];

// 						bool cond =  distance_less_than(iterator1,iterator2,cut_off);
		
// 						if(cond) {
// 							for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							//int q = ints.get_potential_number(iterator1,iterator2,l);
// 							index1_private.push_back(iterator1);
// 							index2_private.push_back(iterator2);
// 							index3_private.push_back(l);	
// 							}				
// 						}						
// 						// for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 						// 	int q = ints.get_potential_number(iterator1,iterator2,l);

// 						// 	if(int_dl[q]) {
// 						// 		bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
// 						// 		if(cond) {
// 						// 			index1_private.push_back(iterator1);
// 						// 			index2_private.push_back(iterator2);
// 						// 			index3_private.push_back(l);					
// 						// 		}

// 						// 	}
// 						// 	else {
// 						// 		index1_private.push_back(iterator1);
// 						// 		index2_private.push_back(iterator2);
// 						// 		index3_private.push_back(l);
// 						// 	}							
// 						// }						
// 					}
// 				}
// 			}
// 			else {

// 			}
// 	}
// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index1.insert(index1.end(),index1_private.begin(),index1_private.end());
// 	}
// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index2.insert(index2.end(),index2_private.begin(),index2_private.end());
// 	}

// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index3.insert(index3.end(),index3_private.begin(),index3_private.end());
// 	}
// 	}
// 	matrix<int> * a = new s_matrix<int>(index1.size(),3);
// 	//s_matrix<int> pairs(index1.size(),3);
// 	for(int i = 0 ; i < (*a).getNsafe() ; i++) {
// 		(*a)(i,0) = index1[i];
// 		(*a)(i,1) = index2[i];
// 		(*a)(i,2) = index3[i];
// 	}

	
	
// 	return a;

// }

// __global__ matrix<int> GPUprecalculatepairs(vector<vector<int> > &b, matrix<int> &boxlist, double cut_off) {
// int n = 0;
// for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++) {
// 	for(int c2 = 0 ; c2 < boxlist.getNsafe() ; c2++) {
// 			int box1 = c1;
// 			int box2 = boxlist(c1,c2);
// 			if(box1==box2) {
// 				n+=(b[box1]).size()*((b[box1]).size()-1)/2;
// 			}
// 			else{
// 				n+=(b[box1].size())*(b[box2].size());
// 			}
// 	}
// }

// int a1(n);
// int a2(n);
// bool c1(n);

// int m = dimension*n;
// double a1pos(m);
// double a2pos(m);

// int iter = 0;
// for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++) {
// 	for(int c2 = 0 ; c2 < boxlist.getNsafe() ; c2++) {
// 		int box1 = c1;
// 		int box2 = boxlist(c1,c2);
// 		if(box1==box2) {
// 			for(int i = 0 ; i < b[box1].size() ; i++ ) {
// 				for(int j = i+1 ;  j < b[box2].size() ; j++) {
// 					int iterator1 = (b[box1])[i];
// 					int iterator2 = (b[box2])[j];
// 					a1[iter]=iterator1;
// 					a2[iter]=iterator2;
// 					iter++;
// 					for(int k = 0  ; k < dimension ; k++) {
// 						a1pos[iter2] = (*dat)(i,k);
// 						a2pos[iter2] = (*dat)(i,k);
// 						iter2++;
// 					}
// 				}
// 			}
// 		}
// 		else{
// 			for(int i = 0 ; i < b[box1].size() ; i++) {
// 				for(int j = 0 ;  j < b[box2].size() ; j++) {
// 					int iterator1 = (b[box1])[i];
// 					int iterator2 = (b[box2])[j];
// 					a1[iter]=iterator1;
// 					a2[iter]=iterator2;
// 					iter++;
// 					for(int k = 0  ; k < dimension ; k++) {
// 						a1pos[iter2] = (*dat)(i,k);
// 						a2pos[iter2] = (*dat)(i,k);
// 						iter2++;
// 					}					
// 		}	
// 	}
// }





// for(int i = 0 ; i < n ; i++) {
// 	c1[i]=false;
// }

// // int iter2 = 0;
// // for(int i = 0 ; i < n ; i++) {
// // 	for(int j = 0  ; j < dimension ; j++) {
// // 		a1pos[i*dimension+j] = (*dat)(i,j);
// // 		a2pos[i*dimension+j] = (*dat)(i,j);
// // 	}
// }

// int *dev_a1,*dev_a2;
// bool *dev_c1;
// double *dev_a1pos,*dev_a2pos;

// HANDLE_ERROR(cudaMalloc((void**)&dev_a1,n*sizeof(int)));
// HANDLE_ERROR(cudaMalloc((void**)&dev_a2,n*sizeof(int)));
// HANDLE_ERROR(cudaMalloc((void**)&dev_c1,n*sizeof(bool)));
// HANDLE_ERROR(cudaMalloc((void**)&dev_a1pos,m*sizeof(double)));
// HANDLE_ERROR(cudaMalloc((void**)&dev_a2pos,m*sizeof(double)));

// HANDLE_ERROR(cudaMemcpy(dev_a1,a1,n*sizeof(int),cudaMemcpyHostToDevice));
// HANDLE_ERROR(cudaMemcpy(dev_a2,a2,n*sizeof(int),cudaMemcpyHostToDevice));
// HANDLE_ERROR(cudaMemcpy(dev_a1pos,a1pos,m*sizeof(double),cudaMemcpyHostToDevice));
// HANDLE_ERROR(cudaMemcpy(dev_a2pos,a2pos,m*sizeof(double),cudaMemcpyHostToDevice));

// //we now have a vector a1,a2 of all the particles which can interact



// }

matrix<int> MD::precalculatepairs(vector<vector<int> > &b, matrix<int> &boxlist, double cut_off) {


	//estimate the total number


	int ss = boxlist.getncols();
	int totn = 0;
	for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++ ) {
		for(int c2  = 0 ; c2 < ss ; c2++) {
			int box1 = c1;
			int box2 = boxlist(c1,c2);
			if(box1==box2) {
				totn += ((b[box1].size()) * (b[box1].size()-1)) / 2;
			}
			else{
				totn += (b[box1]).size() * (b[box2]).size();
			}
		}
	}

	//estimate the total number
	vector<int> index1;
	vector<int> index2;

	index1.reserve(totn);
	index2.reserve(totn);

	//cout << "reserved: " << totn << endl;

	int totb = boxlist.getNsafe();
	#pragma omp parallel
	{
	vector<int> index1_private;
	vector<int> index2_private;

	index1_private.reserve(totn);
	index2_private.reserve(totn);
	//vector<int> index3_private;
	#pragma omp for nowait schedule(static)
	for(int c1 = 0 ; c1 < totb ; c1++)
		for(int c2 = 0 ; c2 < ss ; c2++) {
			int box1 = c1;
			int box2 = boxlist(c1,c2);
			int b1s = b[box1].size();
			int b2s = b[box2].size();

			if(box1==box2) {
				for(int i = 0 ; i < b1s ; i++) {
					for(int j = i+1 ;  j < b2s ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];
						//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);

							//if(int_dl[q]) {
						bool cond = geo.distance_less_than(*dat,iterator1,iterator2,cut_off);
						//bool cond =  distance_less_than(iterator1,iterator2,cut_off);
		
						if(cond) {
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);
							index1_private.push_back(iterator1);
							index2_private.push_back(iterator2);
							//index3_private.push_back(0);	
							//}				
						}							
					}
				}
			}
			else if(box2>box1) {
				for(int i = 0 ; i < b1s ; i++) {
					for(int j = 0 ;  j < b2s ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];

						bool cond = geo.distance_less_than(*dat,iterator1,iterator2,cut_off);
		
						if(cond) {
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);
							index1_private.push_back(iterator1);
							index2_private.push_back(iterator2);
							//index3_private.push_back(0);	
						//	}				
						}						
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
	#pragma omp for schedule(static) ordered
	for(int i = 0 ; i < omp_get_num_threads(); i++) {
		#pragma omp ordered
		index2.insert(index2.end(),index2_private.begin(),index2_private.end());
	}

	}

	matrix<int> a;//(index1.size(),2);
	a.resize_parallel(index1.size(),2);
	//s_matrix<int> pairs(index1.size(),3);
	int NN = (a).getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < NN ; i++) {
		(a)(i,0) = index1[i];
		(a)(i,1) = index2[i];
	}

	
	
	return a;
}

matrix<int> MD::precalculatepairs(const matrix<int> &b, const vector1<int>& sizes, matrix<int> &boxlist, double cut_off)
{

	//estimate the total number

	int ss = boxlist.getncols();
	// int totn = 0;
	// for (int c1 = 0; c1 < boxlist.getNsafe(); c1++)
	// {
	// 	for (int c2 = 0; c2 < ss; c2++)
	// 	{
	// 		int box1 = c1;
	// 		int box2 = boxlist(c1, c2);
	// 		if (box1 == box2)
	// 		{
	// 			totn += ((sizes.gpcons(box1)) * (sizes.gpcons(box1) - 1)) / 2;
	// 		}
	// 		else
	// 		{
	// 			totn += (sizes.gpcons(box1)) * (sizes.gpcons(box2));
	// 		}
	// 	}
	// }
	int N = this->getN();

	int totn = 20*N;


	//estimate the total number
	vector<int> index1;
	vector<int> index2;

	index1.reserve(totn);
	index2.reserve(totn);
	int totb = boxlist.getNsafe();
	//cout << "reserved: " << totn << endl;

#pragma omp parallel
	{
		vector<int> index1_private;
		vector<int> index2_private;

		index1_private.reserve(totn);
		index2_private.reserve(totn);
//vector<int> index3_private;
		#pragma omp for nowait schedule(dynamic)
		for (int c1 = 0; c1 < totb; c1++)
			for (int c2 = 0; c2 < ss; c2++)
			{
				int box1 = c1;
				int box2 = boxlist(c1, c2);
				int b1s = sizes.gpcons(box1);//[box1];
				int b2s = sizes.gpcons(box2);

				if (box1 == box2)
				{
					for (int i = 0; i < b1s; i++)
					{
						for (int j = i + 1; j < b2s; j++)
						{
							int iterator1 = b.gpcons(box1, i);
							int iterator2 = b.gpcons(box2, j);
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);

							//if(int_dl[q]) {
							bool cond = geo.distance_less_than(*dat, iterator1, iterator2, cut_off);
							//bool cond =  distance_less_than(iterator1,iterator2,cut_off);

							if (cond)
							{
								//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
								//int q = ints.get_potential_number(iterator1,iterator2,l);
								index1_private.push_back(iterator1);
								index2_private.push_back(iterator2);
								//index3_private.push_back(0);
								//}
							}
						}
					}
				}
				else if (box2 > box1)
				{
					for (int i = 0; i < b1s; i++)
					{
						for (int j = 0; j < b2s; j++)
						{
							int iterator1 = b.gpcons(box1, i);
							int iterator2 = b.gpcons(box2, j);

							bool cond = geo.distance_less_than(*dat, iterator1, iterator2, cut_off);

							if (cond)
							{
								//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
								//int q = ints.get_potential_number(iterator1,iterator2,l);
								index1_private.push_back(iterator1);
								index2_private.push_back(iterator2);
								//index3_private.push_back(0);
								//	}
							}
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
				else
				{
				}
			}
#pragma omp for schedule(static) ordered
		for (int i = 0; i < omp_get_num_threads(); i++)
		{
#pragma omp ordered
			index1.insert(index1.end(), index1_private.begin(), index1_private.end());
		}
#pragma omp for schedule(static) ordered
		for (int i = 0; i < omp_get_num_threads(); i++)
		{
#pragma omp ordered
			index2.insert(index2.end(), index2_private.begin(), index2_private.end());
		}
	}
	matrix<int> a; //(index1.size(),2);
	a.resize_parallel(index1.size(), 2);
//s_matrix<int> pairs(index1.size(),3);
	int NN = a.getNsafe();


#pragma omp parallel for
	for (int i = 0; i < NN; i++)
	{
		(a)(i, 0) = index1[i];
		(a)(i, 1) = index2[i];
	}

	return a;
}

vector<mdpair> MD::precalculatepairs_md(const matrix<int> &b, const vector1<int> &sizes, matrix<int> &boxlist, double cut_off)
{

	//estimate the total number

	int ss = boxlist.getncols();
	int totn = 0;
	for (int c1 = 0; c1 < boxlist.getNsafe(); c1++)
	{
		for (int c2 = 0; c2 < ss; c2++)
		{
			int box1 = c1;
			int box2 = boxlist(c1, c2);
			if (box1 == box2)
			{
				totn += ((sizes.gpcons(box1)) * (sizes.gpcons(box1) - 1)) / 2;
			}
			else
			{
				totn += (sizes.gpcons(box1)) * (sizes.gpcons(box2));
			}
		}
	}

	//estimate the total number
	vector<mdpair> index1;


	index1.reserve(totn);
	
	int totb = boxlist.getNsafe();
	//cout << "reserved: " << totn << endl;

#pragma omp parallel
	{
		vector<mdpair> index1_private;
		index1_private.reserve(totn);

//vector<int> index3_private;
		#pragma omp for nowait schedule(static)
		for (int c1 = 0; c1 < totb; c1++)
			for (int c2 = 0; c2 < ss; c2++)
			{
				int box1 = c1;
				int box2 = boxlist(c1, c2);
				int b1s = sizes.gpcons(box1); //[box1];
				int b2s = sizes.gpcons(box2);

				if (box1 == box2)
				{
					for (int i = 0; i < b1s; i++)
					{
						for (int j = i + 1; j < b2s; j++)
						{
							int iterator1 = b.gpcons(box1, i);
							int iterator2 = b.gpcons(box2, j);
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);

							//if(int_dl[q]) {
							bool cond = geo.distance_less_than(*dat, iterator1, iterator2, cut_off);
							//bool cond =  distance_less_than(iterator1,iterator2,cut_off);

							if (cond)
							{
								//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
								//int q = ints.get_potential_number(iterator1,iterator2,l);
								index1_private.push_back(mdpair(iterator1,iterator2));
								//index2_private.push_back(iterator2);
								//index3_private.push_back(0);
								//}
							}
						}
					}
				}
				else if (box2 > box1)
				{
					for (int i = 0; i < b1s; i++)
					{
						for (int j = 0; j < b2s; j++)
						{
							int iterator1 = b.gpcons(box1, i);
							int iterator2 = b.gpcons(box2, j);

							bool cond = geo.distance_less_than(*dat, iterator1, iterator2, cut_off);

							if (cond)
							{
								//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
								//int q = ints.get_potential_number(iterator1,iterator2,l);
								index1_private.push_back(mdpair(iterator1, iterator2)); 
								//index2_private.push_back(iterator2);
								//index3_private.push_back(0);
								//	}
							}
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
				else
				{
				}
			}
#pragma omp for schedule(static) ordered
		for (int i = 0; i < omp_get_num_threads(); i++)
		{
#pragma omp ordered
			index1.insert(index1.end(), index1_private.begin(), index1_private.end());
		}

	}
	return index1;
}

matrix<int>* MD::calculatepairs(matrix<int> &boxlist, double cut_off) {
	// vector<int> index1;
	// vector<int> index2;
	//#pragma omp parallel for
	// int Ns = (*dat).getNsafe();
	// int np = ints.number_of_potentials();
	// //vector1<double> int_dis(np);
	// vector1<double> int_dis(np);
	// vector1<bool> int_dl(np);
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
	// }
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dl[i] = (ints.access_potential(i).dl);
	// }


	int total_cubes = boxlist.getNsafe();

	double dims = (double)this->getdimension();

	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

	vector<vector<int> > b;
	int ressize = pow((int)ceil(cut_off),dimension);
	b.reserve(total_cubes);
    for(int j = 0 ; j < total_cubes ; j++) {
        vector<int> temp;
		temp.reserve(ressize);
        b.push_back(temp);
    }
    //int dimension = this->getdimension();


	vector1<int> dim(dimension);
	for(int i = 0 ; i < dimension ; i++ ) {
		int ij = 1;
		for(int j = 0 ; j < i ; j++ ) {
		ij*= cubes_per_length;
		}
		dim[i] = ij;
	}



	for(int i = 0 ; i < this->getN() ; i++) {

		int c = geo.assign_box((*dat),i,dim,cubes_per_length);

		b[c].push_back(i);
	}




	int ss = boxlist.getncols();
	matrix<int> *a = new matrix<int>();
	*a = precalculatepairs(b,boxlist,cut_off);
	return a;


}

matrix<int> *MD::calculatepairs_parallel(matrix<int> &boxlist, double cut_off)
{
	// vector<int> index1;
	// vector<int> index2;
	//#pragma omp parallel for
	// int Ns = (*dat).getNsafe();
	// int np = ints.number_of_potentials();
	// //vector1<double> int_dis(np);
	// vector1<double> int_dis(np);
	// vector1<bool> int_dl(np);

	// for(int i = 0 ; i < np ; i++) {
	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
	// }

	// for(int i = 0 ; i < np ; i++) {
	// 	int_dl[i] = (ints.access_potential(i).dl);
	// }

	int total_cubes = boxlist.getNsafe();

	double dims = (double)this->getdimension();

	int cubes_per_length = (int)round(exp(log(total_cubes) / dims));

	//vector<vector<int>> b;
	int ressize = pow((int)ceil(cut_off), dimension);
	// b.reserve(total_cubes);
	// for (int j = 0; j < total_cubes; j++)
	// {
	// 	vector<int> temp;
	// 	temp.reserve(ressize);
	// 	b.push_back(temp);
	// }
	// //int dimension = this->getdimension();
	// cout << ressize << endl;
	// cout << total_cubes << endl;
	// pausel();

	vector1<int> ccs(total_cubes);
	
	matrix<int> b;

	vector1<int> dim(dimension);
	for (int i = 0; i < dimension; i++)
	{
		int ij = 1;
		for (int j = 0; j < i; j++)
		{
			ij *= cubes_per_length;
		}
		dim[i] = ij;
	}
	

	int partn = this->getN();
	vector1<int> indexes(partn);
	#pragma omp parallel for
	for (int i = 0; i < partn; i++)
	{
		
		int c = geo.assign_box((*dat), i, dim, cubes_per_length);
		indexes[i] = c;
	}



	SingleHistogramParallel(indexes,ccs,b);



	int ss = boxlist.getncols();
	matrix<int> *a = new matrix<int>();
	*a = precalculatepairs(b,ccs, boxlist, cut_off);

	return a;
}


template <class vec>
matrix<int>* MD::calculatepairs(matrix<int> &boxlist, vec &p1, double cut_off) { //p1 is a subset, which interacts with itself
	//ASSUME THE SET OF INTERSECTIONS BETWEEN P1 and P2 is of size 0

	// vector<int> index1;
	// vector<int> index2;
	//#pragma omp parallel for
	// int Ns = (*dat).getNsafe();
	// int np = ints.number_of_potentials();
	// //vector1<double> int_dis(np);
	// vector1<double> int_dis(np);
	// vector1<bool> int_dl(np);
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
	// }
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dl[i] = (ints.access_potential(i).dl);
	// }


	int total_cubes = boxlist.getNsafe();

	double dims = (double)this->getdimension();

	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

	vector<vector<int> > b;
	b.reserve(total_cubes);
    for(int j = 0 ; j < total_cubes ; j++) {
        vector<int> temp;
        b.push_back(temp);
    }
    //int dimension = this->getdimension();


	vector1<int> dim(dimension);
	for(int i = 0 ; i < dimension ; i++ ) {
		int ij = 1;
		for(int j = 0 ; j < i ; j++ ) {
		ij*= cubes_per_length;
		}
		dim[i] = ij;
	}

	// cout << "nonp: " << p1.size() << endl;

	for(int i = 0 ; i < p1.size() ; i++) {

		int c = geo.assign_box((*dat),p1[i],dim,cubes_per_length);

		b[c].push_back(p1[i]);
	}

	//EVEYRTHING BEYOND HERE IS THE PARALLEL ALGO
	
	//pausel();
	//PARALLEL ALGO UP TO HERE

	//int ss = boxlist.getncols();
	matrix<int> *a = new matrix<int>();
	*a = precalculatepairs(b,boxlist,cut_off);
		
	//cout << "nonp: " << a->getnrows() << endl;
	return a;
	

}

template <class vec>
matrix<int> *MD::calculatepairs_parallel(matrix<int> &boxlist, vec &p1, double cut_off)
{ //p1 is a subset, which interacts with itself
	//ASSUME THE SET OF INTERSECTIONS BETWEEN P1 and P2 is of size 0

	// vector<int> index1;
	// vector<int> index2;
	//#pragma omp parallel for
	// int Ns = (*dat).getNsafe();
	// int np = ints.number_of_potentials();
	// //vector1<double> int_dis(np);
	// vector1<double> int_dis(np);
	// vector1<bool> int_dl(np);

	// for(int i = 0 ; i < np ; i++) {
	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
	// }

	// for(int i = 0 ; i < np ; i++) {
	// 	int_dl[i] = (ints.access_potential(i).dl);
	// }

	int total_cubes = boxlist.getNsafe();

	double dims = (double)this->getdimension();

	int cubes_per_length = (int)round(exp(log(total_cubes) / dims));

	int ressize = pow((int)ceil(cut_off), dimension);
	// b.reserve(total_cubes);
	// for (int j = 0; j < total_cubes; j++)
	// {
	// 	vector<int> temp;
	// 	temp.reserve(ressize);
	// 	b.push_back(temp);
	// }
	// //int dimension = this->getdimension();
	// cout << ressize << endl;
	// cout << total_cubes << endl;
	// pausel();

	vector1<int> ccs(total_cubes);

	matrix<int> b;
	//int dimension = this->getdimension();

	vector1<int> dim(dimension);
	for (int i = 0; i < dimension; i++)
	{
		int ij = 1;
		for (int j = 0; j < i; j++)
		{
			ij *= cubes_per_length;
		}
		dim[i] = ij;
	}

	int partn =  p1.size();
	vector1<int> indexes(partn);

	#pragma omp parallel for
	for (int i = 0; i < partn; i++)
	{

		int c = geo.assign_box((*dat), p1[i], dim, cubes_per_length);

		//b[c].push_back(p1[i]);
		indexes[i] = c;
	}



	//when they are added out of order, we do not correspond the entries properly?

	//int ss = boxlist.getncols();
	SingleHistogramParallel(indexes, ccs, b);

	//remap onto the actual values of b
	#pragma omp parallel for
	for(int i = 0 ; i < b.getnrows() ; i++) {
		for(int j = 0  ; j < ccs[i] ; j++) {
			b(i,j) = p1[b(i,j)];
		}
	}


	//find all the particles 

	
	// for (int j = 0; j < b.getnrows(); j++)
	// {
	// 	if (b(j,0) != 0)
	// 	{
	// 		cout << j << " " << endl;
	// 		for (int k = 0; k < b.getncols(); k++)
	// 			cout << b(j,k) << ",";
	// 		cout << endl;
	// 		cout << endl;
	// 	}
	// }

	// pausel();

	int ss = boxlist.getncols();
	matrix<int> *a = new matrix<int>();
	*a = precalculatepairs(b, ccs, boxlist, cut_off);
	


	

	return a;
}

template <class vec>
matrix<int>* MD::calculatepairs(matrix<int> &boxlist, vec &p1, vec &p2, double cut_off) { //p1 and p2 interact with each other but not themselves
	//ASSUME THE SET OF INTERSECTIONS BETWEEN P1 and P2 is of size 0

	vector<int> index1;
	vector<int> index2;
	//#pragma omp parallel for
	// int Ns = (*dat).getNsafe();
	// int np = ints.number_of_potentials();
	// //vector1<double> int_dis(np);
	// vector1<double> int_dis(np);
	// vector1<bool> int_dl(np);
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
	// }
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dl[i] = (ints.access_potential(i).dl);
	// }


	int total_cubes = boxlist.getNsafe();

	double dims = (double)this->getdimension();

	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

	vector<vector<int> > b1;
	vector<vector<int> > b2;
	b1.reserve(total_cubes);
	b2.reserve(total_cubes);
    for(int j = 0 ; j < total_cubes ; j++) {
        vector<int> temp;
        b1.push_back(temp);
        b2.push_back(temp);
    }
    //int dimension = this->getdimension();


	vector1<int> dim(dimension);
	for(int i = 0 ; i < dimension ; i++ ) {
		int ij = 1;
		for(int j = 0 ; j < i ; j++ ) {
		ij*= cubes_per_length;
		}
		dim[i] = ij;
	}



	for(int i = 0 ; i < p1.size() ; i++) {

		int c = geo.assign_box((*dat),p1[i],dim,cubes_per_length);

		b1[c].push_back(p1[i]);
	}
	for(int i = 0 ; i < p2.size() ; i++) {

		int c = geo.assign_box((*dat),p2[i],dim,cubes_per_length);

		b2[c].push_back(p2[i]);
	}

	int tb = boxlist.getNsafe();
	int ss = boxlist.getncols();

	#pragma omp parallel
	{
	vector<int> index1_private;
	vector<int> index2_private;
	//vector<int> index3_private;
	#pragma omp for nowait schedule(static)
	for(int c1 = 0 ; c1 < tb ; c1++) //loop over boxes
		for(int c2 = 0 ; c2 < ss ; c2++) {
			int box1 = c1;
			int box2 = boxlist(c1,c2);
			
			for(int i = 0 ; i < b1[box1].size() ; i++) {
				for(int j = 0 ;  j < b2[box2].size() ; j++) {
					int iterator1 = (b1[box1])[i];
					int iterator2 = (b2[box2])[j];
					//cout << iterator1 << ", " << iterator2 << endl;
					//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
						//int q = ints.get_potential_number(iterator1,iterator2,l);

						//if(int_dl[q]) {
					bool cond = geo.distance_less_than(*dat,iterator1,iterator2,cut_off);
					//bool cond =  distance_less_than(iterator1,iterator2,cut_off);
	
					if(cond) {
						//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
						//int q = ints.get_potential_number(iterator1,iterator2,l);
						index1_private.push_back(iterator1);
						index2_private.push_back(iterator2);
						//index3_private.push_back(0);	
						//}				
					}							
				}
			}
	}
	#pragma omp for schedule(static) ordered
	for(int i = 0 ; i < omp_get_num_threads(); i++) {
		#pragma omp ordered
		index1.insert(index1.end(),index1_private.begin(),index1_private.end());
	}
	#pragma omp for schedule(static) ordered
	for(int i = 0 ; i < omp_get_num_threads(); i++) {
		#pragma omp ordered
		index2.insert(index2.end(),index2_private.begin(),index2_private.end());
	}

	}
	matrix<int> * a = new matrix<int>(index1.size(),2);
	//s_matrix<int> pairs(index1.size(),3);
	for(int i = 0 ; i < (*a).getNsafe() ; i++) {
		(*a)(i,0) = index1[i];
		(*a)(i,1) = index2[i];
	}

	
	
	return a;


}


matrix<int>* MD::calculatepairs_sorted(matrix<int> &boxlist, double cut_off) {
	vector<mdpair> index1;
	//vector<int> index2;
	//#pragma omp parallel for
	// int Ns = (*dat).getNsafe();
	// int np = ints.number_of_potentials();
	// //vector1<double> int_dis(np);
	// vector1<double> int_dis(np);
	// vector1<bool> int_dl(np);
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
	// }
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dl[i] = (ints.access_potential(i).dl);
	// }


	int total_cubes = boxlist.getNsafe();

	double dims = (double)this->getdimension();

	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

	vector<vector<int> > b;
	b.reserve(total_cubes);
    for(int j = 0 ; j < total_cubes ; j++) {
        vector<int> temp;
        b.push_back(temp);
    }
    //int dimension = this->getdimension();


	vector1<int> dim(dimension);
	for(int i = 0 ; i < dimension ; i++ ) {
		int ij = 1;
		for(int j = 0 ; j < i ; j++ ) {
		ij*= cubes_per_length;
		}
		dim[i] = ij;
	}



	for(int i = 0 ; i < this->getN() ; i++) {

		int c = geo.assign_box((*dat),i,dim,cubes_per_length);

		b[c].push_back(i);
	}


	int tb =boxlist.getNsafe();
	int ss = boxlist.getncols();


	#pragma omp parallel
	{
	vector<mdpair> index1_private;
	//vector<int> index2_private;
	//vector<int> index3_private;
	#pragma omp for nowait schedule(static)
	for(int c1 = 0 ; c1 < tb ; c1++)
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
						bool cond = geo.distance_less_than(*dat,iterator1,iterator2,cut_off);
						//bool cond =  distance_less_than(iterator1,iterator2,cut_off);
		
						if(cond) {
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);
							index1_private.push_back(mdpair(iterator1,iterator2));
							//index2_private.push_back(iterator2);
							//index3_private.push_back(0);	
							//}				
						}							
					}
				}
			}
			else if(box2>box1) {
				for(int i = 0 ; i < b[box1].size() ; i++) {
					for(int j = 0 ;  j < b[box2].size() ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];

						bool cond = geo.distance_less_than(*dat,iterator1,iterator2,cut_off);
		
						if(cond) {
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);
							index1_private.push_back(mdpair(iterator1,iterator2));
							//index3_private.push_back(0);	
						//	}				
						}						
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
	// #pragma omp for schedule(static) ordered
	// for(int i = 0 ; i < omp_get_num_threads(); i++) {
	// 	#pragma omp ordered
	// 	index2.insert(index2.end(),index2_private.begin(),index2_private.end());
	// }
	}
	sort(index1.begin(),index1.end());
	matrix<int> * a = new matrix<int>(index1.size(),2);
	//s_matrix<int> pairs(index1.size(),3);
	for(int i = 0 ; i < (*a).getNsafe() ; i++) {
		(*a)(i,0) = index1[i].a;
		(*a)(i,1) = index1[i].b;
	}

	
	
	return a;


}

matrix<double> MD::calculateforces(matrix<int> &pairs) {
	
	matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	int totp = pairs.getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < totp ; i++ ) {
		int p1 = pairs.mat[i*2+0];
		int p2 = pairs.mat[i*2+1];
		//int i1 = pairs(i,2);

		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		geo.distance_vector(*dat,p1,p2,un,dis);
		double f1  = (*ints).force(sqrt(dis));

		for(int j = 0 ; j < dimension ; j++) {
			double fac = f1*un[j]/sqrt(dis);
			(forces).mat[p1*dimension+j] += fac;
			(forces).mat[p2*dimension+j] += -fac;
		}
	}



	return forces;
}

matrix<double> MD::calculateforces(matrix<int> &pairs,potential &iny) {
	
	matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();



	int totp = pairs.getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < totp ; ++i ) {
		int p1 = pairs.mat[i*2+0];
		int p2 = pairs.mat[i*2+1];
		//int i1 = pairs(i,2);
		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		geo.distance_vector(*dat,p1,p2,un,dis);

		//un = i-j


		double f1  = iny.force(sqrt(dis));
	
		if(abs(f1) > 1.E4) {
			cout << "force too big" << endl;

			cout << p1 << " " << p2 << endl;
			cout << f1 << " " << dis << endl;
			vector1<int> dim(3);
			for (int i1 = 0; i1 < 3; i1++)
			{
				int ij = 1;
				for (int j = 0; j < i1; j++)
				{
					ij *= 38;
				}
				dim[i] = ij;
			}

			vector1<double> pr1 = (*dat).getrowvector(p1);

			vector1<double> pr2 = (*dat).getrowvector(p2);

			vector1<int> pri1(3);
			vector1<int> pri2(3);

			for (size_t i2 = 0; i2 < 3; i2++)
			{
				pri1[i2] = floor(pr1[i2] / 3.05908);
				pri2[i2] = floor(pr2[i2] / 3.05908);
				
			}
			

			cout << pr1 << endl;
			cout << pr2 << endl;

			cout << scalar(dim,pri1) << endl;
			cout << scalar(dim,pri2) << endl;
			
			pausel();
		}  
		

		

		for(int j = 0 ; j < dimension ; ++j) {
			double fac = f1*un[j]/sqrt(dis);
			(forces).mat[p1*dimension+j] += fac;
			(forces).mat[p2*dimension+j] += -fac;
		}
	}



	return forces;
}

matrix<double> MD::calculateforces(vector<mdpair> &pairs,potential &iny) {
//NOT the most efficient way to do this, but I am being quick right now

matrix<int> pp(pairs.size(),2);
for(int i = 0  ; i < pairs.size() ; i++) {
	pp(i, 0) = pairs[i].a;
	pp(i, 1) = pairs[i].b;
}
return calculateforces(pp,iny);

}

matrix<double> MD::calculateforcesDV(matrix<int> &pairs, potential &iny, vector1<double> &mags)
{

	matrix<double> forces((*dat).getNsafe(), dimension);
	// vec_vec<double> outputs(pairs.getn());
	//  ofstream myfile;
	//  myfile.open("forces.csv", ios::out | ios::app);
	// cout << pairs.getNsafe() << endl;
	// potential * pot = ints(0,1,0).clone();
	int totp = pairs.getNsafe();
#pragma omp parallel for
	for (int i = 0; i < totp; ++i)
	{
		int p1 = pairs.mat[i * 2 + 0];
		int p2 = pairs.mat[i * 2 + 1];
		// int i1 = pairs(i,2);
		double dis;
		// vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		geo.distance_vector(*dat, p1, p2, un, dis);

		// un = i-j

		double f1 = mags[i]*iny.force(sqrt(dis));
		/*
			if(abs(f1) > 1.E4) {
				cout << p1 << " " << p2 << endl;
				cout << f1 << " " << dis << endl;
				vector1<int> dim(3);
				for (int i1 = 0; i1 < 3; i1++)
				{
					int ij = 1;
					for (int j = 0; j < i1; j++)
					{
						ij *= 38;
					}
					dim[i] = ij;
				}

				vector1<double> pr1 = (*dat).getrowvector(p1);

				vector1<double> pr2 = (*dat).getrowvector(p2);

				vector1<int> pri1(3);
				vector1<int> pri2(3);

				for (size_t i2 = 0; i2 < 3; i2++)
				{
					pri1[i2] = floor(pr1[i2] / 3.05908);
					pri2[i2] = floor(pr2[i2] / 3.05908);

				}


				cout << pr1 << endl;
				cout << pr2 << endl;

				cout << scalar(dim,pri1) << endl;
				cout << scalar(dim,pri2) << endl;

				pausel();
			} */

		for (int j = 0; j < dimension; ++j)
		{
			double fac = f1 * un[j] / sqrt(dis);
			(forces).mat[p1 * dimension + j] += fac;
			(forces).mat[p2 * dimension + j] += -fac;
		}
	}

	return forces;
}

matrix<double> MD::calculateforces_sp(matrix<int> &pairs, potential &iny, matrix<vector1<double> > &sp)
{

	matrix<double> forces((*dat).getNsafe(), dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	int totp = pairs.getNsafe();
#pragma omp parallel for
	for (int i = 0; i < totp; ++i)
	{
		int p1 = pairs.mat[i * 2 + 0];
		int p2 = pairs.mat[i * 2 + 1];
		//int i1 = pairs(i,2);
		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		geo.distance_vector(*dat, p1, p2, un, dis);

		//un = i-j
		iny.setparameters(sp(p1,p2));
		double f1 = iny.force(sqrt(dis));

		for (int j = 0; j < dimension; ++j)
		{
			double fac = f1 * un[j] / sqrt(dis);
			(forces).mat[p1 * dimension + j] += fac;
			(forces).mat[p2 * dimension + j] += -fac;
		}
	}

	return forces;
}


matrix<double> MD::calculateforceslist(matrix<int> &pairs,potential &iny) {
	
	matrix<double> forces(pairs.getNsafe(),2*dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	int totp = pairs.getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < totp ; ++i ) {
		int p1 = pairs.mat[i*2+0];
		int p2 = pairs.mat[i*2+1];
		//int i1 = pairs(i,2);
		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		geo.distance_vector(*dat,p1,p2,un,dis);

		//un = i-j

		double f1  = iny.force(sqrt(dis));

		for(int j = 0 ; j < dimension ; ++j) {
			double fac = f1*un[j]/sqrt(dis);
			// (forces).mat[p1*dimension+j] += fac;
			// (forces).mat[p2*dimension+j] += -fac;
			forces(i,j)=fac;
			forces(i,j+dimension)=-fac;
		}
	}



	return forces;
}




matrix<double> MD::calculateforces_threebody(matrix<int> &triplets,potential3 &iny) {
	
	matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	int totp = triplets.getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < totp ; ++i ) {
		int p1 = triplets.mat[i*3+0];
		int p2 = triplets.mat[i*3+1];
		int p3 = triplets.mat[i*3+2];
		//int i1 = pairs(i,2);
		double dis1,dis2;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un1(dimension),un2(dimension);
		geo.distance_vector(*dat,p2,p1,un1,dis1);
		geo.distance_vector(*dat,p3,p2,un2,dis2);

		vector1<double> f1(dimension);
		vector1<double> f2(dimension);
		vector1<double> f3(dimension);
		iny.force(un1,un2,f1,f2,f3);

		// if(chckmatrix(fas)) {
		// 	cout << un1 << endl;
		// 	cout << un2 << endl;
		// 	cout << p1 << endl;
		// 	cout << p2 << endl;
		// 	cout << p3 << endl;
		// 	cout << fas << endl;
		// 	error("error in bending forces");
		// }

		for(int k = 0 ; k < dimension ; k++) {
		forces(p1,k) += f1[k];
		forces(p2,k) += f2[k];
		forces(p3,k) += f3[k];
		}

		// double f1  = iny.force(sqrt(dis));

		// for(int j = 0 ; j < dimension ; j++) {
		// 	double fac = f1*un[j]/sqrt(dis);
		// 	(forces).mat[p1*dimension+j] += fac;
		// 	(forces).mat[p2*dimension+j] += -fac;
		// }
	}



	return forces;
}

matrix<double> MD::calculateforces_threebody(vector<mdtriplet> &trips, potential3 &iny) {
	matrix<int> pp(trips.size(),3);
	for (int i = 0; i < trips.size(); i++)
	{
		pp(i, 0) = trips[i].a;
		pp(i, 1) = trips[i].b;
		pp(i, 2) = trips[i].c;

		
	}
	return calculateforces_threebody(pp, iny);
}

matrix<double> MD::calculateforces_fast3D(matrix<int> &pairs) {
	
	matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	int totp = pairs.getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < totp ; i++ ) {


		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		//int p1 = pairs.mat[i*2+0];
		//int p2 = pairs.mat[i*2+1];
		//geo.distance_vector(*dat,p1,p2,un,dis);

		int p1= pairs.mat[i*2+0]*dimension;
		int p2 =pairs.mat[i*2+1]*dimension;

		geo.distance3v((*dat).mat[p1+0],(*dat).mat[p2+0],(*dat).mat[p1+1],(*dat).mat[p2+1],(*dat).mat[p1+2],(*dat).mat[p2+2],un,dis);

		double f1  = (*ints).force2mdx(dis);

		for(int j = 0 ; j < dimension ; j++) {
			double fac = f1*un[j];
			(forces).mat[pairs.mat[i*2+0]*dimension+j] += fac;
			(forces).mat[pairs.mat[i*2+1]*dimension+j] += -fac;
		}
	}



	return forces;
}

matrix<double> MD::calculateforcesharmonic(matrix<int> &pairs, matrix<double> &bondlengths, double k)
{

	matrix<double> forces((*dat).getNsafe(), dimension);
	// vec_vec<double> outputs(pairs.getn());
	//  ofstream myfile;
	//  myfile.open("forces.csv", ios::out | ios::app);
	// cout << pairs.getNsafe() << endl;
	// potential * pot = ints(0,1,0).clone();

	int totp = pairs.getNsafe();
#pragma omp parallel for
	for (int i = 0; i < totp; ++i)
	{
		int p1 = pairs.mat[i * 2 + 0];
		int p2 = pairs.mat[i * 2 + 1];
		// int i1 = pairs(i,2);
		double dis;
		// vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		geo.distance_vector(*dat, p1, p2, un, dis);

		// un = i-j

		double f1 = -k*(sqrt(dis)-bondlengths[i]);//iny.force(sqrt(dis));

		if (abs(f1) > 1.E4)
		{
			cout << "force too big" << endl;

			cout << p1 << " " << p2 << endl;
			cout << f1 << " " << dis << endl;
			vector1<int> dim(3);
			for (int i1 = 0; i1 < 3; i1++)
			{
				int ij = 1;
				for (int j = 0; j < i1; j++)
				{
					ij *= 38;
				}
				dim[i] = ij;
			}

			vector1<double> pr1 = (*dat).getrowvector(p1);

			vector1<double> pr2 = (*dat).getrowvector(p2);

			vector1<int> pri1(3);
			vector1<int> pri2(3);

			for (size_t i2 = 0; i2 < 3; i2++)
			{
				pri1[i2] = floor(pr1[i2] / 3.05908);
				pri2[i2] = floor(pr2[i2] / 3.05908);
			}

			cout << pr1 << endl;
			cout << pr2 << endl;

			cout << scalar(dim, pri1) << endl;
			cout << scalar(dim, pri2) << endl;

			pausel();
		}

		for (int j = 0; j < dimension; ++j)
		{
			double fac = f1 * un[j] / sqrt(dis);
			(forces).mat[p1 * dimension + j] += fac;
			(forces).mat[p2 * dimension + j] += -fac;
		}
	}

	return forces;
}

matrix<double> MD::calculateforcesdelauny(matrix<int> &quads, double kappa) {

	matrix<double> forces((*dat).getNsafe(), dimension);
	int totp = quads.getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < totp ; i++ ) {
		int p1 = quads.mat[i * 4 + 0];
		int p2 = quads.mat[i * 4 + 1];
		int p3 = quads.mat[i * 4 + 2];
		int p4 = quads.mat[i * 4 + 3];

		double x1 = (*dat)[p1 * dimension + 0];
		double y1 = (*dat)[p1 * dimension + 1];
		double z1 = (*dat)[p1 * dimension + 2];

		double x2 = (*dat)[p2 * dimension + 0];
		double y2 = (*dat)[p2 * dimension + 1];
		double z2 = (*dat)[p2 * dimension + 2];

		double x3 = (*dat)[p3 * dimension + 0];
		double y3 = (*dat)[p3 * dimension + 1];
		double z3 = (*dat)[p3 * dimension + 2];

		double x4 = (*dat)[p4 * dimension + 0];
		double y4 = (*dat)[p4 * dimension + 1];
		double z4 = (*dat)[p4 * dimension + 2];

		// double x1 = 0.;
		// double y1 = 0.;
		// double z1 = 0.;

		// double x2 = 1.;
		// double y2 = 0.;
		// double z2 = 0.;

		// double x3 = 0.5;
		// double y3 = sqrt(3)/2;
		// double z3 = -.3;

		// double x4 = 0.5;
		// double y4 = -sqrt(3)/2;;
		// double z4 = -.3;

		//assume that very sharp local fluctuations are impossible

		//both normal vectors should be pointing out of the plane of the vesicle, or into the plane of the veiscle
		double g1 = SQR(-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + SQR((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) + SQR(-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3)) ;
		double g2 = SQR((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + SQR(-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4)) + SQR((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4)) ;
		double normn = sqrt(g1) * sqrt(g2);
		double fac1= sqrt(g2)/sqrt(g1);
		double fac2 = 1./fac1;
		double sumn = (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) +
																																																					   ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4)) +
																																																					   (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3)) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));

		// double dnormdx1;
		// double dnormdy1;
		// double dnormdz1;
		// double dnormdx2;
		// double dnormdy2;
		// double dnormdz2;
		// double dnormdx3;
		// double dnormdy3;
		// double dnormdz3;
		// double dnormdx4;
		// double dnormdy4;
		// double dnormdz4;
		double dg1dx1 = 2 * (y2 - y3) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + 2 * (-z2 + z3) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3));

		double dg2dx1 = 2 * (-y2 + y4) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + 2 * (z2 - z4) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4));

		double dg1dy1 = 2 * (-x2 + x3) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + 2 * (z2 - z3) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3));

		double dg2dy1 = 2 * (x2 - x4) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + 2 * (-z2 + z4) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));

		double dg1dz1 = 2 * (x2 - x3) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) + 2 * (-y2 + y3) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3));

		double dg2dz1 = 2 * (-x2 + x4) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4)) + 2 * (y2 - y4) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));

		double dg1dx2 = 2 * (-y1 + y3) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + 2 * (z1 - z3) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3));

		double dg2dx2 = 2 * (y1 - y4) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + 2 * (-z1 + z4) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4));

		double dg1dy2 = 2 * (x1 - x3) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + 2 * (-z1 + z3) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3));

		double dg2dy2 = 2 * (-x1 + x4) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + 2 * (z1 - z4) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));

		double dg1dz2 = 2 * (-x1 + x3) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) + 2 * (y1 - y3) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3));

		double dg2dz2 = 2 * (x1 - x4) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4)) + 2 * (-y1 + y4) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));

		double dg1dx3 = 2 * (y1 - y2) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + 2 * (-z1 + z2) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3));

		double dg2dx3 = 0;

		double dg1dy3 = 2 * (-x1 + x2) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + 2 * (z1 - z2) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3));

		double dg2dy3 = 0;

		double dg1dz3 = 2 * (x1 - x2) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) + 2 * (-y1 + y2) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3));

		double dg2dz3 = 0;

		double dg1dx4 = 0;

		double dg2dx4 = 2 * (-y1 + y2) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + 2 * (z1 - z2) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4));

		double dg1dy4 = 0;

		double dg2dy4 = 2 * (x1 - x2) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + 2 * (-z1 + z2) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));

		double dg1dz4 = 0;

		double dg2dz4 = 2 * (-x1 + x2) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4)) + 2 * (y1 - y2) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));


		double dnormdx1 = 0.5 * (fac1 * dg1dx1 + fac2 * dg2dx1);
		double dnormdy1 = 0.5 * (fac1 * dg1dy1 + fac2 * dg2dy1);
		double dnormdz1 = 0.5 * (fac1 * dg1dz1 + fac2 * dg2dz1);
		double dnormdx2 = 0.5 * (fac1 * dg1dx2 + fac2 * dg2dx2);
		double dnormdy2 = 0.5 * (fac1 * dg1dy2 + fac2 * dg2dy2);
		double dnormdz2 = 0.5 * (fac1 * dg1dz2 + fac2 * dg2dz2);
		double dnormdx3 = 0.5 * (fac1 * dg1dx3 + fac2 * dg2dx3);
		double dnormdy3 = 0.5 * (fac1 * dg1dy3 + fac2 * dg2dy3);
		double dnormdz3 = 0.5 * (fac1 * dg1dz3 + fac2 * dg2dz3);
		double dnormdx4 = 0.5 * (fac1 * dg1dx4 + fac2 * dg2dx4);
		double dnormdy4 = 0.5 * (fac1 * dg1dy4 + fac2 * dg2dy4);
		double dnormdz4 = 0.5 * (fac1 * dg1dz4 + fac2 * dg2dz4);

		double dsumndx1 = (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) * (-y2 + y4) + (y2 - y3) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) * (z2 - z4) + (-z2 + z3) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4));
		double dsumndy1 = (x2 - x4) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + (-x2 + x3) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3)) * (-z2 + z4) + (z2 - z3) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));
		double dsumndz1 = (-x2 + x4) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) + (y2 - y4) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3)) + (x2 - x3) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4)) + (-y2 + y3) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));
		double dsumndx2 = (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) * (y1 - y4) + (-y1 + y3) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) * (-z1 + z4) + (z1 - z3) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4));
		double dsumndy2 = (-x1 + x4) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + (x1 - x3) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3)) * (z1 - z4) + (-z1 + z3) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));
		double dsumndz2 = (x1 - x4) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) + (-y1 + y4) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3)) + (-x1 + x3) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4)) + (y1 - y3) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));
		double dsumndx3 = (y1 - y2) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + (-z1 + z2) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4));
		double dsumndy3 = (-x1 + x2) * ((-x2 + x4) * (-y1 + y4) - (-x1 + x4) * (-y2 + y4)) + (z1 - z2) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));
		double dsumndz3 = (x1 - x2) * (-((-x2 + x4) * (-z1 + z4)) + (-x1 + x4) * (-z2 + z4)) + (-y1 + y2) * ((-y2 + y4) * (-z1 + z4) - (-y1 + y4) * (-z2 + z4));
		double dsumndx4 = (-y1 + y2) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + (z1 - z2) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3));
		double dsumndy4 = (x1 - x2) * (-((-x1 + x3) * (-y1 + y2)) + (-x1 + x2) * (-y1 + y3)) + (-z1 + z2) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3));
		double dsumndz4 = (-x1 + x2) * ((-x1 + x3) * (-z1 + z2) - (-x1 + x2) * (-z1 + z3)) + (y1 - y2) * (-((-y1 + y3) * (-z1 + z2)) + (-y1 + y2) * (-z1 + z3));

		double f1 = kappa * (sumn * dnormdx1 - normn * dsumndx1) / (2 * SQR(normn));
		double f2 = kappa *(sumn * dnormdy1 - normn * dsumndy1) / (2 * SQR(normn));
		double f3 = kappa *(sumn * dnormdz1 - normn * dsumndz1) / (2 * SQR(normn));

		double f4 = kappa * (sumn * dnormdx2 - normn * dsumndx2) / (2 * SQR(normn));
		double f5 = kappa * (sumn * dnormdy2 - normn * dsumndy2) / (2 * SQR(normn));
		double f6 = kappa * (sumn * dnormdz2 - normn * dsumndz2) / (2 * SQR(normn));

		double f7 = kappa * (sumn * dnormdx3 - normn * dsumndx3) / (2 * SQR(normn));
		double f8 = kappa * (sumn * dnormdy3 - normn * dsumndy3) / (2 * SQR(normn));
		double f9 = kappa * (sumn * dnormdz3 - normn * dsumndz3) / (2 * SQR(normn));

		double f10 = kappa * (sumn * dnormdx4 - normn * dsumndx4) / (2 * SQR(normn));
		double f11 = kappa * (sumn * dnormdy4 - normn * dsumndy4) / (2 * SQR(normn));
		double f12 = kappa * (sumn * dnormdz4 - normn * dsumndz4) / (2 * SQR(normn));

		// bool cond1 =  abs(f1) > 1E3 || abs(f2) > 1E3 || abs(f3) > 1E3 || abs(f4) > 1E3 || abs(f5) > 1E3 || abs(f6) > 1E3 || abs(f7) > 1E3 || abs(f8) > 1E3 || abs(f9) > 1E3 || abs(f10) > 1E3 || abs(f11) > 1E3 || abs(f12) > 1E3;
		// if(cond1) {
		// 	cout << x1 << " " << y1 << " " << z1 << endl;
		// 	cout << x2 << " " << y2 << " " << z2 << endl;
		// 	cout << x3 << " " << y3 << " " << z3 << endl;
		// 	cout << x4 << " " << y4 << " " << z4 << endl;
		// 	cout << f1 << endl;
		// 	cout << f2 << endl;
		// 	cout << f3 << endl;
		// 	cout << f4 << endl;
		// 	cout << f5 << endl;
		// 	cout << f6 << endl;
		// 	cout << f7 << endl;
		// 	cout << f8 << endl;
		// 	cout << f9 << endl;
		// 	cout << f10 << endl;
		// 	cout << f11 << endl;
		// 	cout << f12 << endl;
		// 	cout << normn << endl;
		// 	cout << sumn << endl;
		// 	pausel();
		// }

		// double v1[3] = {(y3 * (z1 - z2) + y1 * (z2 - z3) + y2 * (-z1 + z3)) / (SQR(x1) + x2 * x3 - x1 * (x2 + x3) + (y1 - y2) * (y1 - y3) + (z1 - z2) * (z1 - z3)),
		// 				(x3 * (-z1 + z2) + x2 * (z1 - z3) + x1 * (-z2 + z3)) / (SQR(x1) + x2 * x3 - x1 * (x2 + x3) + (y1 - y2) * (y1 - y3) + (z1 - z2) * (z1 - z3)),
		// 				(x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3)) / (SQR(x1) + x2 * x3 - x1 * (x2 + x3) + (y1 - y2) * (y1 - y3) + (z1 - z2) * (z1 - z3)) };

		// double v2[3] = {(y4 * (-z1 + z2) + y2 * (z1 - z4) + y1 * (-z2 + z4)) / ((-x1 + x4) * (-x2 + x4) + (-y1 + y4) * (-y2 + y4) + (-z1 + z4) * (-z2 + z4)), (x4 * (z1 - z2) + x1 * (z2 - z4) + x2 * (-z1 + z4)) / ((-x1 + x4) * (-x2 + x4) + (-y1 + y4) * (-y2 + y4) + (-z1 + z4) * (-z2 + z4)), (x4 * (-y1 + y2) + x2 * (y1 - y4) + x1 * (-y2 + y4)) / ((-x1 + x4) * (-x2 + x4) + (-y1 + y4) * (-y2 + y4) + (-z1 + z4) * (-z2 + z4))};

		// double v3[3] = {(x1 + x2 + x3)/3.,(y1 + y2 + y3)/3.,(z1 + z2 + z3)/3.};

		// double v4[3] = {(x1 + x2 + x4) / 3., (y1 + y2 + y4) / 3., (z1 + z2 + z4) / 3.};

		// matrix<double> m(12,3);
		// m(0, 0) = v1[0];
		// m(0, 1) = v1[1];
		// m(0, 2) = v1[2];

		// m(1, 0) = v2[0];
		// m(1, 1) = v2[1];
		// m(1, 2) = v2[2];

		// m(2, 0) = v3[0];
		// m(2, 1) = v3[1];
		// m(2, 2) = v3[2];

		// m(3, 0) = v4[0];
		// m(3, 1) = v4[1];
		// m(3, 2) = v4[2];

		// m(4, 0) = x1;
		// m(4, 1) = y1;
		// m(4, 2) = z1;

		// m(5, 0) = x2;
		// m(5, 1) = y2;
		// m(5, 2) = z2;

		// m(6, 0) = x3;
		// m(6, 1) = y3;
		// m(6, 2) = z3;

		// m(7, 0) = x4;
		// m(7, 1) = y4;
		// m(7, 2) = z4;

		// m(8, 0) = f1;
		// m(8, 1) = f2;
		// m(8, 2) = f3;

		// m(9, 0) = f4;
		// m(9, 1) = f5;
		// m(9, 2) = f6;

		// m(10, 0) = f7;
		// m(10, 1) = f8;
		// m(10, 2) = f9;

		// m(11, 0) = f10;
		// m(11, 1) = f11;
		// m(11, 2) = f12;

		// outfunc(m,"test");
		// pausel();



		forces.mat[p1 * dimension + 0] += -f1;
		forces.mat[p1 * dimension + 1] += -f2;
		forces.mat[p1 * dimension + 2] += -f3;

		forces.mat[p2 * dimension + 0] += -f4;
		forces.mat[p2 * dimension + 1] += -f5;
		forces.mat[p2 * dimension + 2] += -f6;

		forces.mat[p3 * dimension + 0] += -f7;
		forces.mat[p3 * dimension + 1] += -f8;
		forces.mat[p3 * dimension + 2] += -f9;

		forces.mat[p4 * dimension + 0] += -f10;
		forces.mat[p4 * dimension + 1] += -f11;
		forces.mat[p4 * dimension + 2] += -f12;


	}


	return forces;
}

matrix<double> MD::calculatestress(matrix<int> &pairs) {
	
	//matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	matrix<double> stress(dimension,dimension);
	int totp = pairs.getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < totp ; i++ ) {


		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		//int p1 = pairs.mat[i*2+0];
		//int p2 = pairs.mat[i*2+1];
		//geo.distance_vector(*dat,p1,p2,un,dis);

		int p1= pairs.mat[i*2+0]*dimension;
		int p2 =pairs.mat[i*2+1]*dimension;

		geo.distance3v((*dat).mat[p1+0],(*dat).mat[p2+0],(*dat).mat[p1+1],(*dat).mat[p2+1],(*dat).mat[p1+2],(*dat).mat[p2+2],un,dis);

		double f1  = (*ints).force2mdx(dis);

		for(int k = 0 ; k < dimension ; k++) {
		for(int j = 0 ; j < dimension ; j++) {
			
			double fac = f1*sqrt(dis)*un[k]*un[j];
			stress(k,j) += fac;
			// double fac2 = 
			// (forces).mat[pairs.mat[i*2+0]*dimension+j] += fac;
			// (forces).mat[pairs.mat[i*2+1]*dimension+j] += -fac;
		}
		}
	}



	return stress;
}

matrix<double> MD::calculateforces_truncated(matrix<int> &pairs, double above) {
	
	matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	int totp = pairs.getNsafe();
	#pragma omp parallel for
	for(int i = 0 ; i < totp ; i++ ) {


		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		//int p1 = pairs.mat[i*2+0];
		//int p2 = pairs.mat[i*2+1];
		//geo.distance_vector(*dat,p1,p2,un,dis);

		int p1= pairs.mat[i*2+0]*dimension;
		int p2 =pairs.mat[i*2+1]*dimension;

		geo.distance3v((*dat).mat[p1+0],(*dat).mat[p2+0],(*dat).mat[p1+1],(*dat).mat[p2+1],(*dat).mat[p1+2],(*dat).mat[p2+2],un,dis);

		double f1  = (*ints).force2mdx(dis);
		if(dis < above) f1 = 0;

		for(int j = 0 ; j < dimension ; j++) {
			double fac = f1*un[j];
			(forces).mat[pairs.mat[i*2+0]*dimension+j] += fac;
			(forces).mat[pairs.mat[i*2+1]*dimension+j] += -fac;
		}
	}



	return forces;
}


class spherical_confinement_3D {
	private:
	double rmax;

	double v;

	double shft;


	public:
	spherical_confinement_3D() : rmax(10.0), v(1.0), shft(5.)
	{
	}

	spherical_confinement_3D(double rmaxx, double vv, double shftt) : rmax(rmaxx), v(vv), shft(shftt)
	{

	}
	spherical_confinement_3D& operator=(const spherical_confinement_3D &old) {
		rmax = old.rmax;
		v = old.v;
		shft = old.shft;
		return *this ;
	}

	void setv(double vv) {
		v=vv;
	}

	vector1<double> operator()(const vector1<double> &x) const {
		
		//for(int j = 0  ; j < 3 ; )
		// double r = sqrt(SQR(x.gpcons(0)) + SQR(x.gpcons(1)) + SQR(x.gpcons(2)));
		double r = sqrt(SQR(x.gpcons(0)-shft) + SQR(x.gpcons(1)-shft) + SQR(x.gpcons(2)-shft));



		vector1<double> force(3);

		if(abs(r-rmax) > 2.) return force;
		else {

		double r2 = SQR(r-rmax);

		double r6 = CUB(r2);

		double r12 = SQR(r6);

		force[0] = 12 * v * (x.gpcons(0)-shft) / (r * (r - rmax) * r12);
		force[1] = 12 * v * (x.gpcons(1)-shft) / (r * (r - rmax) * r12);
		force[2] = 12 * v * (x.gpcons(2)-shft) / (r * (r - rmax) * r12);

		return force;
		}
	}

};

class planar_confinement {
	private :
		vector1<double> lmax;
		double v;

	public:
		planar_confinement() : lmax(vector1<double>(3,10.)), v(1.0)
		{
		}

		planar_confinement(vector1<double> rmaxx, double vv) : lmax(rmaxx), v(vv)
		{
		}
		planar_confinement &operator=(const planar_confinement &old)
		{
			lmax = old.lmax;
			v = old.v;
			return *this;
		}
		void setv(double vv)
		{
			v = vv;
		}

		void setl(const vector1<double> &ll) {
			lmax = ll;
		}

		vector1<double> operator()(const vector1<double> &x) const
		{

			// for(int j = 0  ; j < 3 ; )
			//  double r = sqrt(SQR(x.gpcons(0)) + SQR(x.gpcons(1)) + SQR(x.gpcons(2)));
			//double r = sqrt(SQR(x.gpcons(0) - shft) + SQR(x.gpcons(1) - shft) + SQR(x.gpcons(2) - shft));

			vector1<double> force(3);


				for(int j = 0 ; j < 3 ; j++) {
				double r2 = SQR(x.gpcons(j) - lmax.gpcons(j));

				double r6 = CUB(r2);

				double r12 = SQR(r6);

				double q2 = SQR(x.gpcons(j));

				double q6 = CUB(r2);

				double q12 = SQR(r6);

				force[j] = v * ((12. / (x.gpcons(j) * q12)) + (12. / ((x.gpcons(2)-lmax.gpcons(j)) * r12)));
				}
				return force;
			
		}
};


template <class Q>
matrix<double> MD::calculateforces_external(Q &func)
{

	matrix<double> forces((*dat).getNsafe(), dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	int totp = (*dat).getNsafe();
	#pragma omp parallel for
	for (int i = 0; i < totp; i++)
	{

		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un = (*dat).getrowvector(i);
		//int p1 = pairs.mat[i*2+0];
		//int p2 = pairs.mat[i*2+1];
		//geo.distance_vector(*dat,p1,p2,un,dis);

		// int p1 = pairs.mat[i * 2 + 0] * dimension;
		// int p2 = pairs.mat[i * 2 + 1] * dimension;


		vector1<double> f =  func(un);
		// cout << f << endl;
		// pausel();
		// geo.distance3v((*dat).mat[p1 + 0], (*dat).mat[p2 + 0], (*dat).mat[p1 + 1], (*dat).mat[p2 + 1], (*dat).mat[p1 + 2], (*dat).mat[p2 + 2], un, dis);

		//double f1 = (*ints).force2mdx(dis);

		for (int j = 0; j < dimension; j++)
		{
			//double fac = f1 * un[j];
			//if(f[j]>0) { cout << f[j] << endl; pausel(); }
			(forces).mat[i* dimension + j] += f[j];
			//(forces).mat[pairs.mat[i * 2 + 1] * dimension + j] += -fac;
		}
	}

	return forces;
}

//void MD::adv(matrix<int> &pairs) { error("adv function not overriden matrix");}