Microtubule::Microtubule(double ll, int Na, int Nb, int nmic, int LL) : pai(vector1<int>(Na)),pbi(vector1<int>(Nb)),pci(vector1<int>(LL*nmic)),bound(vector1<int>(Na+Nb)),bound_along(vector1<double>(Na+Nb)) {
	obj = new LangevinNVT;
	dimension = 2;
	l = ll;
	is_periodic = true;
	num = floor(ll/4.0);
	double pi = 3.1415;

	// chm = 0.5*l;
	spatial_dependence =  false;
	double sigma = 1.0;
	double distance_between_points = 2.0;
	L = LL;
	Ld = (double)L;
	vector1<bool> pb(dimension,true);
	cube bc(ll,pb,dimension);

	na = Na;
	nb = Nb;
	nc = nmic*L;
	number_of_microtubules = nmic;

	for(int i = 0 ; i < na ; i++)
		pai[i]=i;
	int iter2 = 0;
	for(int i = na ; i < na+nb ; i++) {
		pbi[iter2]=i;
		iter2++;
	}

	int iter3 = 0;
	for(int i = na+nb ; i < na+nb+nc ; i++) {
		pci[iter3]=i;
		iter3++;
	}

	for(int i = 0 ; i < na+nb ; i++)  {
		bound[i]=0; //no binding to begin with
		bound_along[i]=0.0;
	}

	totalN = Na+Nb+nc;

	


	matrix<double> res(nmic,dimension-1); // orientation vector of the microtubules
	for(int i = 0 ; i < number_of_microtubules ; i++) {
		res(i,0)=pi/2.;//2*pi*(double)rand()/(double)RAND_MAX;
	}


	orient  = new matrix<double>;
	*orient = res;
//	intmatrix f(totalN);
	
	WCAPotential nhj(2.0,sigma,2.0);
	potential* q1 = nhj.clone();
	faa = q1;

	potential* q2 = nhj.clone();
	fbb = q2;


	WCAPotential nswca(10.0, sigma, 0.0);
	potential* q3 = nswca.clone();
	fcc = q3;

	potential* q4 = nswca.clone();
	fab = q4;

	potential* q5 = nswca.clone();
	fbc = q5;

	potential* q6 = nswca.clone();
	fac = q6;	

	HarmonicPotential nfr(100.0,1.0);
	potential* q7 = nfr.clone();
	bindp = q7;	


	FENEPotential nrr(50.0,1.2);
	potential *q8 = nrr.clone();
	bindm =  q8;

	BendingPotential nfr2(1000.0,0.0);
	potential3* q9 = nfr2.clone();
	bendp = q9;
	//double specific_volume  = sigma*sigma*sigma;
	//double dls = sigma;

	int pp = floor(l/sigma);
	vector<double> possible_pos_x;
	vector<double> possible_pos_y; 
	
	
	// possible_pos_x_rods.reserve(pp*pp);
	// possible_pos_y_rods.reserve(pp*pp);

	possible_pos_x.reserve(pp*pp);
	possible_pos_y.reserve(pp*pp);
	
	for(int i = 0 ; i < pp ; i++)
		for(int j = 0 ; j < pp ; j++) {
				double x =  0.5*sigma+i*sigma;
				double y =  0.5*sigma+j*sigma;
				//vector<double> b;
				possible_pos_x.push_back(x);
				possible_pos_y.push_back(y);
				//possible_pos.push_back(b);
	}
		
	//pausel();
		

	cout << "system possible set up complete" << endl;

	// double xh = 4.0;
	// double yh = 50.0;
	// bool ghj1 = xh>L/2. ? true : false ;
	// bool ghj2 = xh<l-L/2. ? true : false ;
	// bool ghj3 = (yh>L/2.) ? true : false ;
	// bool ghj4 = yh<l-L/2. ? true : false ;
 // 	cout << ghj1 << endl;
 // 	cout << ghj2 << endl;
 // 	cout << ghj3 << endl;
 // 	cout << ghj4 << endl;
 // 	cout << (ghj1&&ghj2&&ghj3&&ghj4) << endl;
 // 	pausel();
	//cout << totalN << endl;
	matrix<double> store(totalN,dimension);
	int gh = (int)((double)(L-1.)/(2.*sigma));
	//double b = 1.0;
	//cout << gh << endl;
//	cout << gh << endl;
	//double orient  = 0;
	for(int i = 0 ; i < number_of_microtubules ; i++) { //SET UP POLYMERS
		//cout << i << endl;
		int randint =(rand() % (possible_pos_x.size())); //position of middle point
		do{
		randint = (rand() % (possible_pos_x.size()));
		//for(int j = 0 ; j < dimension ; j++)
		for(int j = 0 ; j <= 2*gh ; j++) {
		
		//cout << j << " " << i*L+na+nb+j << endl;
		store(i*L+na+nb+j,0) =  possible_pos_x[randint];
		store(i*L+na+nb+j,1) =  possible_pos_y[randint]-gh+j;
		}
		// cout << i << endl;
		// cout << possible_pos_x[randint] << ", " << possible_pos_y[randint] << endl; 
		//cout << randint << endl;
		//cout << possible_pos_x[randint] << ", " << possible_pos_y[randint] << endl;
		//store(i,1) =  possible_pos[randint][1];
		//store(i,2) =  possible_pos[randint][2];
		// for(int j = 0 ; j < 2*gh+1 ; j++) {
		// //cout << j << endl;
		// 	possible_pos_x


		} while ( !(possible_pos_x[randint]>L? true : false) && !(possible_pos_x[randint]<l-L ? true : false) && !(possible_pos_y[randint]>L ? true : false) && !(possible_pos_y[randint]<l-L ? true : false) );
		// }
		//cout << possible_pos_x.size() << endl;

		// cout << randint << endl;
		// for(int k = randint-gh ;  k <= randint+gh ; k++) {
		// cout << possible_pos_x[k] << " " << possible_pos_y[k] << endl;
		// }
		// for(int k =na+nb ; k < totalN ; k++) {
		// cout << store(k,'r') << endl;
		// }
		// pausel();

		vector<double>::iterator it1,it2;
		it1 = possible_pos_x.begin();
		it1 += randint - gh;
		it2 = possible_pos_x.begin();
		it2 += randint + gh+1;

		vector<double>::iterator it3,it4;
		it3 = possible_pos_y.begin();
		it3 += randint - gh;
		it4 = possible_pos_y.begin();
		it4 += randint + gh+1;		
		possible_pos_x.erase(it1,it2);
		possible_pos_y.erase(it3,it4);

	}	


		// for(int j = 0 ; j < dimension ; j++) {
		// store(i*L+na+nb,0) =  possible_pos_x[randint];
		// store(i*L+na+nb,1) =  possible_pos_y[randint];
		// }
		// for(int j = 1 ; j < L ; j++) {
		// 	store(i*L+na+nb+j,0) = possible_pos_x[randint]+sigma*cos(orient);
		// 	store(i*L+na+nb+j,1) = possible_pos_y[randint]+sigma*sin(orient);
		// }
	

	// for(int i = 0 ; i < nc ; i++) { SET UP FOR LINES
	// 	//cout << i << endl;
	// 	int randint =(rand() % (possible_pos_x.size()));
	// 	do{
	// 	randint = (rand() % (possible_pos_x.size()));
	// 	//for(int j = 0 ; j < dimension ; j++)
	// 	store(i+na+nb,0) =  possible_pos_x[randint];
	// 	store(i+na+nb,1) =  possible_pos_y[randint];
	// 	// cout << i << endl;
	// 	// cout << possible_pos_x[randint] << ", " << possible_pos_y[randint] << endl; 
	// 	//cout << randint << endl;
	// 	//cout << possible_pos_x[randint] << ", " << possible_pos_y[randint] << endl;
	// 	//store(i,1) =  possible_pos[randint][1];
	// 	//store(i,2) =  possible_pos[randint][2];
	// 	// for(int j = 0 ; j < 2*gh+1 ; j++) {
	// 	// //cout << j << endl;
	// 	// 	possible_pos_x


	// 	} while ( !(possible_pos_x[randint]>L/2.? true : false) && !(possible_pos_x[randint]<l-L/2. ? true : false) && !(possible_pos_y[randint]>L/2. ? true : false) && !(possible_pos_y[randint]<l-L/2. ? true : false) );
	// 	// }
	// 	//cout << possible_pos_x.size() << endl;
	// 	vector<double>::iterator it1,it2;
	// 	it1 = possible_pos_x.begin();
	// 	it1 += randint - gh;
	// 	it2 = possible_pos_x.begin();
	// 	it2 += randint + gh;

	// 	vector<double>::iterator it3,it4;
	// 	it3 = possible_pos_y.begin();
	// 	it3 += randint - gh;
	// 	it4 = possible_pos_y.begin();
	// 	it4 += randint + gh;		
	// 	possible_pos_x.erase(it1,it2);
	// 	possible_pos_y.erase(it3,it4);

	// }

	// cout << possible_pos_x.size() << endl;
	//pausel();


	for(int i = 0 ; i < na+nb ; i++) {
		//cout << i << endl;
		int randint = rand() % (possible_pos_x.size());
		//for(int j = 0 ; j < dimension ; j++)
		store(i,0) =  possible_pos_x[randint];
		store(i,1) =  possible_pos_y[randint];
		//store(i,1) =  possible_pos[randint][1];
		//store(i,2) =  possible_pos[randint][2];
		possible_pos_x.erase(possible_pos_x.begin()+randint);
		possible_pos_y.erase(possible_pos_y.begin()+randint);
	}	

cout << "system actual set up complete" << endl;



matrix<int> bondpairss(nmic*(L-1),2);

int iter5 = 0;
for(int i = 0  ; i < number_of_microtubules ; i++) {
	for(int j = 0 ; j < L-1 ; j++) {
		//cout << iter5 << endl;
		bondpairss(iter5,0)=i*L+na+nb+j;
		bondpairss(iter5,1)=i*L+na+nb+j+1;
		iter5++;
	}
}



bondpairs  = new matrix<int>;
*bondpairs = bondpairss;

matrix<int> bendpairss(nmic*(L-2),3);
int iter6 = 0;
for(int i = 0  ; i < number_of_microtubules ; i++) {
	for(int j = 0 ; j < L-2 ; j++) {
		//cout  << iter6 << endl;
		bendpairss(iter6,0)=i*L+na+nb+j;
		bendpairss(iter6,1)=i*L+na+nb+j+1;
		bendpairss(iter6,2)=i*L+na+nb+j+2;
		iter6++;
	}
}

cout << 2 << endl;

bendtriplets  = new matrix<int>;
*bendtriplets = bendpairss;
	//cout << store << endl;
// potential *faa;
// potential *fbb;
// potential *fcc;
// potential *fab;
// potential *fbc;
// potential *fac;

	//LJPotential nswca(1.0,1.0);
	//ExponentialPotential nswca(eqeps,1.0);
	
//	f.add_interaction_all(nswca);

	LangevinNVT b(bc);

	double kT = 1.0;
	dt = 0.005;
	double eta = 50.;
	gamma  = eta;

	probunbind = dt*0.00;
	probbind = 1.0;
	
	disless = 1.2;
	excess_force_distance = 1.2;
	v0_a = 33.0;
	v0_b = 33.0;

	polarity_a = 1.0;
	polarity_b = 1.0;

	max_s=SQR(2.);


	//chckcollisions(store);

	
	b.setdat(store);
	b.setinteractions(nswca);
	b.setkT(kT);
	b.setdt(dt);
	b.setgamma(eta);
	b.setm(1.0);
	matrix<double> moms(b.getN(),dimension);
	b.setmom(moms);


	*obj  = b;

	// possible_pos_x.clear();

	// possible_pos_y.clear();

	cout << "created" << endl;
//	pausel();
	// cout << *bondpairs << endl;
	// cout << *bendtriplets << endl;
	// pausel();

}

Microtubule::~Microtubule() {
	delete obj;
	delete bondpairs;
	delete bendtriplets;
	delete orient;
	delete faa;
	delete fbb;
	delete fcc;
	delete fab;
	delete fbc;
	delete fac;
	delete bindp;
	delete bindm;
	delete bendp;
}

void Microtubule::setkT(double kT) {
	obj->setkT(kT);
}

void Microtubule::setgeo(cube &bc) {
	(*(this->obj)).setgeometry(bc);
	l=bc.l;
	num = floor(bc.l/4.0);
}

void Microtubule::setpotaa(potential &a) {
	potential* q = a.clone();
	faa = q;
}

void Microtubule::setpotab(potential &a) {
	potential* q = a.clone();
	fab = q;
}

void Microtubule::setpotac(potential &a) {
	potential* q = a.clone();
	fac = q;
}

void Microtubule::setpotbb(potential &a) {
	potential* q = a.clone();
	fbb = q;
}

void Microtubule::setpotbc(potential &a) {
	potential* q = a.clone();
	fbc = q;
}

void Microtubule::setpotcc(potential &a) {
	potential* q = a.clone();
	fcc = q;
}

void Microtubule::setprobunbind(double a) {
	probunbind = a;
}

void Microtubule::setprobbind(double a) {
	probbind = a;
}

void Microtubule::set_excess_force_distance(double a) {
	excess_force_distance = a;
}

void Microtubule::setv0(double a, double b) {
	v0_a = a;
	v0_b =b;
}

void Microtubule::setpolarity(double a, double b) {
	polarity_a = a;
	polarity_b =b;
}

void Microtubule::setgamma(double a) {
	gamma = a;
	obj->setgamma(gamma);
}

void Microtubule::set_gamma_spatial_dependence(bool temp) {
	spatial_dependence = temp;
}

void Microtubule::set_initial_conditions(string filename1, string filename2, string filename3) {
	double T;
	bool err1=false;
	bool err2=false;
	bool err3=false;
	matrix<double> particleA(na,dimension);
	matrix<double> particleB(nb,dimension);
	matrix<double> microt(nc,dimension);

	if(na>0) {
	particleA = importcsv(filename1,T,err1); //import particles A
	}
	if(nb>0) {
	particleB = importcsv(filename2,T,err2); //import particles B
	}
	if(nc > 0) {
	microt = importcsv(filename3,T,err3); //import microtubules
	}

	if(particleA.getNsafe() != na ) error("size of data in filename A not correct");
	if(particleB.getNsafe() != nb ) error("size of data in filename B not correct");
	if(microt.getNsafe() != nc ) error("size of data in filename C not correct");

	if(err1) error("error in importing file 1");
	if(err2) error("error in importing file 2");
	if(err3) error("error in importing file 3");

	matrix<double> store(totalN,dimension);
	for(int i =  0 ; i < na ; i++) {
		for(int j = 0 ; j < dimension ; j++) {
			store(i,j) = particleA(i,j);
		}
	}
	for(int i =  0 ; i < nb ; i++) {
		for(int j = 0 ; j < dimension ; j++) {
			store(i+na,j) = particleB(i,j);
		}
	}
	for(int i =  0 ; i < nc ; i++) {
		for(int j = 0 ; j < dimension ; j++) {
			store(i+na+nb,j) = microt(i,j);
		}
	}			

	(*obj).setdat(store);
}


matrix<double> Microtubule::PositionForcesDueToAngles() {
	int totalN = obj->getN();
	//what kind of position forces are generated due to angular evolution.
	matrix<double> forces(totalN,dimension);
	//one example is motility:
	vector<double> motorsbound(number_of_microtubules,0.);

	// vector1<double> avuva(dimension);

	// vector1<double> avuvb(dimension);

	for(int i = 0 ; i < na+nb ; i++) {
		if(bound[i]>0) {
			motorsbound[bound[i]-1]+=1.;
			//int p2 = i;
			double tt  =  bound_along[i];
		//	cout << tt << endl;
			int mt =  bound[i];
			double tt2 = tt*(double)(L-1);
		//	cout << tt2 << endl;

			//int k1 = floor(tt2);
			int k2 = ceil(tt2);
			int k1 = k2-1;

			//cout << k1 << endl;
		//	cout << k2 << endl;

			//double t = tt2-k1;

			//cout << t << endl;


			//na+nb+mt
			int p1a = na+nb+(mt-1)*L+k1;
			int p1b = na+nb+(mt-1)*L+k2;

			// cout << p1a << endl;
			// cout << p1b << endl;

			vector1<double> a1(dimension);
			vector1<double> a2(dimension);
			// a0[0]=obj->getcoordinate(p2,0);
			// a0[1]=obj->getcoordinate(p2,1);
			a1[0]=obj->getcoordinate(p1a,0);
			a1[1]=obj->getcoordinate(p1a,1);
			a2[0]=obj->getcoordinate(p1b,0);
			a2[1]=obj->getcoordinate(p1b,1);

			vector1<double> uv(dimension);
			double d;
			(obj->getgeo()).distance_vector(a1, a2, uv, d);

			// cout << uv << endl;
			// cout << d << endl;
			uv /= sqrt(d);

			// cout << bound_along[i] << endl;
			//double incr = SQR(dt)*(1./(1+0.5*gamma*dt))*v0;

			if(i < na) {
				//avuva += polarity_a*uv;
				double incr = dt*v0_a/gamma;
				bound_along[i] = bound_along[i] - polarity_a*incr/((double)(L-1));

				if(bound_along[i]<0) bound_along[i]=0;

				for(int j = 0 ; j < L ; j++) {
					 forces(na+nb+(mt-1)*L+j,0)+=-polarity_a*0.5*v0_a*uv[0]; 
					 forces(na+nb+(mt-1)*L+j,1)+=-polarity_a*0.5*v0_a*uv[1];
				}
			}
			else {
				//avuvb += polarity_b*uv;
				double incr = dt*v0_b/gamma;
				bound_along[i] = bound_along[i] - polarity_b*incr/((double)(L-1));

				if(bound_along[i]<0) bound_along[i]=0;

				for(int j = 0 ; j < L ; j++) {
					 forces(na+nb+(mt-1)*L+j,0)+=-polarity_b*0.5*v0_b*uv[0]; 
					 forces(na+nb+(mt-1)*L+j,1)+=-polarity_b*0.5*v0_b*uv[1];
				}
			}			
			// forces(p2,0)=0.5*v0*uv[0];
			// forces(p2,1)=0.5*v0*uv[1];
			//F(i,2) = v0*z;
			//cout << a1 << endl;
			//cout << a2 << endl;
	//		cout << v0 << endl;
			//cout << uv << endl;
			// cout << forces << endl;
			// pausel();

		}
	}

	// cout << avuva+avuvb << endl;
	// cout << "prenorm: " << forces(na+nb,0) << " " << forces(na+nb,1) << endl;

	//cout << forces(na+nb,0) << "," << forces(na+nb,1) << endl;


	double v0 = (v0_a+v0_b)/2.;
	for(int i = 0 ; i < number_of_microtubules ; i++) { //normalize for many motors (not fully collective)
		if(motorsbound[i]>1E-10) {
			for(int j = 0 ; j < L ; j++) {
				double mag =0.0;
				for(int k = 0 ; k < dimension ; k++) {
					mag+=SQR(forces(na+nb+i*L+j,k));
				}
//				cout << "mag: " << mag << endl;
				double fac = (mag/(max_s*v0*v0))/tanh(mag/(max_s*v0*v0));
//				cout << "fac: " << fac << endl;
				if(fac > 1E-10) {
				double rescale = 1./fac;
				// cout << "rescaled: " << rescale << endl;
				// pausel();
				for(int k = 0 ; k < dimension ; k++) {
					forces(na+nb+i*L+j,k)=rescale*forces(na+nb+i*L+j,k);
				}
				}
			}
		}
	}

	//cout << forces(na+nb,0) << "," << forces(na+nb,1) << endl;
//	cout << "postnorm: " << forces(na+nb,0) << " " << forces(na+nb,1) << endl;

	//cout << forces << endl;
	// for(int i = 0 ; i < nc ; i++) { //this would choose a random new orientation as a perturbation of the previous one
	// 	double theta  = (*orient)(i,0);
	// 	//double phi = (*orient)(i,1);	
	// 	double x = cos(theta);
	// 	double y = sin(theta);
	// 	//double z = cos(theta);

	// 	F(i,0) = v0*x;
	// 	F(i,1) = v0*y;
	// 	//F(i,2) = v0*z;
	// }
	return forces;


}

matrix<double> Microtubule::constantMTforce() //This function constantly propels every Microtubule
{
	int totalN = obj->getN();
	//what kind of position forces are generated due to angular evolution.
	matrix<double> forces(totalN, dimension);
	//one example is motility:

	// vector1<double> avuva(dimension);

	// vector1<double> avuvb(dimension);

	for (int i = 0; i < number_of_microtubules; i++)
	{
		vector1<double> UV(dimension);
		for(int j = 0 ; j < L-1 ; j++)
		{
			vector1<double> uv(dimension);
			double d;
			vector1<double> a1(dimension);
			vector1<double> a2(dimension);

			int p1a = j;
			int p1b = j+1;
			a1[0] = obj->getcoordinate(p1a, 0);
			a1[1] = obj->getcoordinate(p1a, 1);
			a2[0] = obj->getcoordinate(p1b, 0);
			a2[1] = obj->getcoordinate(p1b, 1);
			(obj->getgeo()).distance_vector(a1, a2, uv, d);
			UV += uv;
			}
			UV /= (double)(L-1);
		
		for (int j = 0; j < L; j++)
		{
			forces(na + nb + (i) * L + j, 0) += -polarity_a * 0.5 * v0_a * UV[0];
			forces(na + nb + (i) * L + j, 1) += -polarity_a * 0.5 * v0_a * UV[1];
		}
	}
		return forces;
}


//add constant force to everything
matrix<double> Microtubule::constantMTforce(vector1<int> &p1, vector1<int> &p2) //This function constantly propels every Microtubule
{
	int totalN = obj->getN();
	//what kind of position forces are generated due to angular evolution.
	matrix<double> forces(totalN, dimension);
	//one example is motility:
	
	// vector1<double> avuva(dimension);

	// vector1<double> avuvb(dimension);
	vector1<double> UV(dimension);

	for (int i = 0; i < p1.getsize(); i++)
	{
		

			vector1<double> uv(dimension);
			double d;
			vector1<double> a1(dimension);
			vector1<double> a2(dimension);

			int p1a = p1[i];
			int p1b = p2[i];
			a1[0] = obj->getcoordinate(p1a, 0);
			a1[1] = obj->getcoordinate(p1a, 1);
			a2[0] = obj->getcoordinate(p1b, 0);
			a2[1] = obj->getcoordinate(p1b, 1);
			(obj->getgeo()).distance_vector(a1, a2, uv, d);
			UV += uv;
		

	}
	UV /= (double)(p1.getsize());
	for (int i = 0; i < number_of_microtubules; i++)
	{
		for (int j = 0; j < L; j++)
		{
			forces(na + nb + (i)*L + j, 0) += -polarity_a * 0.5 * v0_a * UV[0];
			forces(na + nb + (i)*L + j, 1) += -polarity_a * 0.5 * v0_a * UV[1];
		}
	}
	return forces;
}

/*
matrix<double> Microtubule::PositionForcesDueToAngles(double v0) {
	int totalN = obj->getN();
	//what kind of position forces are generated due to angular evolution.
	matrix<double> forces(totalN,dimension);
	//one example is motility:
	for(int i = 0 ; i < na+nb ; i++) {
		if(bound[i]>0) {
			int p2 = i;
			double tt  =  bound_along[i];
		//	cout << tt << endl;
			int mt =  bound[i];
			double tt2 = tt*(double)(L-1);
		//	cout << tt2 << endl;

			//int k1 = floor(tt2);
			int k2 = ceil(tt2);
			int k1 = k2-1;

			//cout << k1 << endl;
		//	cout << k2 << endl;

			double t = tt2-k1;

			//cout << t << endl;

			//na+nb+mt
			int p1a = na+nb+(mt-1)*L+k1;
			int p1b = na+nb+(mt-1)*L+k2;

			// cout << p1a << endl;
			// cout << p1b << endl;

			vector1<double> a1(dimension);
			vector1<double> a2(dimension);
			// a0[0]=obj->getcoordinate(p2,0);
			// a0[1]=obj->getcoordinate(p2,1);
			a1[0]=obj->getcoordinate(p1a,0);
			a1[1]=obj->getcoordinate(p1a,1);
			a2[0]=obj->getcoordinate(p1b,0);
			a2[1]=obj->getcoordinate(p1b,1);

			vector1<double> uv(dimension);
			double d;
			(obj->getgeo()).distance_vector(a1, a2, uv, d);

			// cout << uv << endl;
			// cout << d << endl;
			uv /= sqrt(d);

			// cout << bound_along[i] << endl;
			//double incr = SQR(dt)*(1./(1+0.5*gamma*dt))*v0;
			double incr = dt*v0/gamma;

			bound_along[i] = bound_along[i] - incr/((double)(L-1));

			if(bound_along[i]<0) bound_along[i]=0;

			//for(int j = 0 ; j < L ; j++) {
				 forces(k2,0)+=-v0*uv[0]; 
				 forces(k2,1)+=-v0*uv[1];
				 forces(k1,0)+=-v0*uv[0]; 
				 forces(k1,1)+=-v0*uv[1];
			//}
			forces(p2,0)=v0*uv[0];
			forces(p2,1)=v0*uv[1];
			//F(i,2) = v0*z;
			//cout << a1 << endl;
			//cout << a2 << endl;
	//		cout << v0 << endl;
			//cout << uv << endl;
			// cout << forces << endl;
			// pausel();

		}
	}
	
	// for(int i = 0 ; i < nc ; i++) { //this would choose a random new orientation as a perturbation of the previous one
	// 	double theta  = (*orient)(i,0);
	// 	//double phi = (*orient)(i,1);	
	// 	double x = cos(theta);
	// 	double y = sin(theta);
	// 	//double z = cos(theta);

	// 	F(i,0) = v0*x;
	// 	F(i,1) = v0*y;
	// 	//F(i,2) = v0*z;
	// }
	return forces;


}
*/

void Microtubule::minimal_distance_between_point_and_line(vector1<double> &point1, vector1<double> &end1, vector1<double> &end2, double &dis, vector1<double> &uv, double &tt) {
		// where t is the distance along the line

		vector1<double> uv1(dimension);
		vector1<double> uv2(dimension);

		double dis1 = 0.0;
		double dis2 = 0.0;

		(obj->getgeo()).distance_vector(end2,end1,uv1,dis1);
		(obj->getgeo()).distance_vector(end1,point1,uv2,dis2);

		double t =  -scalar(uv1,uv2)/dis1;

		vector1<double> uv3(dimension);
		double dis3=0.0;

		if(t<0||t>1) { //endpoint condition 
			(obj->getgeo()).distance_vector(end2,point1,uv3,dis3);
			if(dis3<dis2) t= 1;
			else t = 0;
		}

		

		vector1<double> a3(dimension);

		a3[0]=end1[0]+t*(uv1[0]);
		a3[1]=end1[1]+t*(uv1[1]);

		double dis4=0.0;
		vector1<double> uv4(dimension);

		(obj->getgeo()).correct_position(a3);

		(obj->getgeo()).distance_vector(a3,point1,uv4,dis4);

		uv = uv4;
		dis =  dis4;
		tt = t-0.5;
}

void Microtubule::minimal_distance_between_point_and_line(int &p1a,int &p1b, int &p2, double &dis, vector1<double> &uv, double &tt) {
		// where t is the distance along the line

		
		vector1<double> point1(dimension);
		vector1<double> end1(dimension);
		vector1<double> end2(dimension);


		//double theta1 = (*orient)(p1-na-nb,0);

		point1[0]=obj->getcoordinate(p2,0);
		point1[1]=obj->getcoordinate(p2,1);

		end1[0]=obj->getcoordinate(p1b,0);//+(L/2.)*cos(theta1);
		end1[1]=obj->getcoordinate(p1b,1);//+(L/2.)*sin(theta1);

		end2[0]=obj->getcoordinate(p1a,0);//-(L/2.)*cos(theta1);
		end2[1]=obj->getcoordinate(p1a,1);//-(L/2.)*sin(theta1);	

		(obj->getgeo()).correct_position(end1);
		(obj->getgeo()).correct_position(end2);

		this->minimal_distance_between_point_and_line(point1, end1, end2, dis, uv, tt);


}


void Microtubule::Bind(int &p1, int &p2, int initialbind, int &finalbind, double &fpos) { //p1 is the microtubule
		int mt = 1+(int)((double)(p1-na-nb)/(double)(L)); //which microtubule, starting from 1
	//	cout << 1 << endl;
		double tt = (double)((p1-na-nb)%L);
	//	cout << 2 << endl;
		tt = tt/(double)(L-1);
	//	cout << 3 << endl;
		double dis4 = obj->distance(p1,p2);
	//	cout << 4 << endl;
		if(initialbind==0 && abs(tt-0.5)<0.4 && dis4 < disless ) { //only if they are close and the particle is close
		double rr1 = ((double) rand() / (RAND_MAX));
			if(rr1<probbind) { //bind with rate
			finalbind=mt; //p2 is bound to p1;
			fpos = tt;
		//	cout << p1 << " " << p2 << " " << tt << endl;
			}
		}
}



void Microtubule::CalculateBindings(matrix<int>&pairs1,matrix<int>&pairs2) { //pairs1 is the distance between type A and pairs2 is the distance between type B
	//bound is a vector of bindings between particles and the microtubule they are attached too.
	//matrix<double> 

	//the same process as calculating the forces.
	vector1<int> tempbound(bound);
	vector1<double> tempbound_along(bound_along);
	// if(trap(bound,1)>0) {
	// 	cout << "before binding calc" << endl;
	// 	cout << bound << endl;
	// 	cout << tempbound << endl;
	// 	cout << bound_along << endl;
	// 	cout << tempbound_along << endl;
	// 	pausel();
	// }
	vector1<bool> unbound(na+nb,false);

		for(int i = 0 ; i < na+nb ; i++) {
			if(bound[i]>0) { 
			//unbinding rate
			int p2 = i;
				//calcualte distance between microtubule:
			double tt  =  bound_along[i];
			int mt =  bound[i];
			double tt2 = tt*(double)(L-1);


			
			int k2 = ceil(tt2);
			int k1 = k2 - 1;
			double t = tt2-k1;
			int p1a = na+nb+(mt-1)*L+k1;
			int p1b = na+nb+(mt-1)*L+k2;
			vector1<double> a0(dimension),a1(dimension),a2(dimension);
			a0[0]=obj->getcoordinate(p2,0);
			a0[1]=obj->getcoordinate(p2,1);
			a1[0]=obj->getcoordinate(p1a,0);
			a1[1]=obj->getcoordinate(p1a,1);
			a2[0]=obj->getcoordinate(p1b,0);
			a2[1]=obj->getcoordinate(p1b,1);
			vector1<double> uv3(dimension),uv4(dimension);
			double dis3,dis4;
			(obj->getgeo()).distance_vector(a2,a1,uv3,dis3);
			vector1<double> tempvector = a1+(t*(uv3));
			(obj->getgeo()).correct_position(tempvector);
			(obj->getgeo()).distance_vector(a0,tempvector,uv4,dis4);

			double rr1 = ((double) rand() / (RAND_MAX));
				if(abs(bound_along[i]-0.5)>0.4) { //drop off
					unbound[i] = true;
					tempbound[i]=0;
					tempbound_along[i]=0.;
				}
				else if(rr1<probunbind) { //unbind with rate
				unbound[i] = true;
				tempbound[i]=0;
				tempbound_along[i]=0.;
				}
				else if(sqrt(dis4)>excess_force_distance) { // if distance is larger
				unbound[i] = true;
				tempbound[i]=0;
				tempbound_along[i]=0.;					
				}
				else {

				}

			}
		}
		//UNBINDINGS 
	// if(trap(bound,1)>0) {
	// 	cout << "after unbinding" << endl;
	// 	cout << bound << endl;
	// 	cout << tempbound << endl;
	// 	cout << bound_along << endl;
	// 	cout << tempbound_along << endl;
	// 	pausel();
	// }
		//cout << "calc" << endl;
		//calculate bindings

//for each pair which could be bound
	for(int i = 0 ; i < pairs1.getNsafe() ; i++ ) {
		int p1 = pairs1(i,0) > pairs1(i,1) ?  pairs1(i,0) : pairs1(i,1) ;
		int p2 = pairs1(i,0) < pairs1(i,1) ?  pairs1(i,0) : pairs1(i,1) ;
		int fb = bound[p2];
		double ba = bound_along[p2];
		if(!unbound[p2]&&tempbound[p2]==0) {
		this->Bind(p1,p2,bound[p2],fb,ba); //only if unbound
		tempbound[p2] =  fb;
		tempbound_along[p2] = ba;
		}

	}

	// if(trap(bound,1)>0) {
	// 	cout << "after bind 1" << endl;
	// 	cout << bound << endl;
	// 	cout << tempbound << endl;
	// 	cout << bound_along << endl;
	// 	cout << tempbound_along << endl;
	// 	pausel();
	// }

	for(int i = 0 ; i < pairs2.getNsafe() ; i++ ) {
		int p1 = pairs2(i,0) > pairs2(i,1) ?  pairs2(i,0) : pairs2(i,1) ;
		int p2 = pairs2(i,0) < pairs2(i,1) ?  pairs2(i,0) : pairs2(i,1) ;
		int fb = bound[p2];
		double ba = bound_along[p2];

		if(!unbound[p2]&&tempbound[p2]==0) {
		this->Bind(p1,p2,bound[p2],fb,ba);		
		tempbound[p2] =  fb;
		tempbound_along[p2] = ba;
		}

	}

	// if(trap(bound,1)>0) {
	// 	cout << "after bind 2" << endl;
	// 	cout << bound << endl;
	// 	cout << tempbound << endl;
	// 	cout << bound_along << endl;
	// 	cout << tempbound_along << endl;
	// 	pausel();
	// }

	bound = tempbound;
	bound_along = tempbound_along;

	//cout << "calc done" << endl;

}

matrix<double> Microtubule::BindingForcesTest(vector1<int> &bound1, vector1<double> &bound_along1, matrix<double> &dat1) {
		matrix<double> forces(totalN,dimension);

		for(int i = 0 ; i < na+nb ; i++ ) {
		int p2 = i;
		int p1 = bound1[i];
		if(p1>0) {

			double tt  =  bound_along1[i];
			int mt =  bound1[i];
			double tt2 = tt*(double)(L-1);


			
			int k2 = ceil(tt2);
			int k1 = k2 - 1;

			double t = tt2-k1;

			//na+nb+mt
			int p1a = na+nb+(mt-1)*L+k1;
			int p1b = na+nb+(mt-1)*L+k2;


			vector1<double> a0(dimension),a1(dimension),a2(dimension);

			// a0[0]=obj->getcoordinate(p2,0);
			// a0[1]=obj->getcoordinate(p2,1);
			// a1[0]=obj->getcoordinate(p1a,0);
			// a1[1]=obj->getcoordinate(p1a,1);
			// a2[0]=obj->getcoordinate(p1b,0);
			// a2[1]=obj->getcoordinate(p1b,1);

			a0[0]=dat1(p2,0);//obj->getcoordinate(p2,0);
			a0[1]=dat1(p2,1);//obj->getcoordinate(p2,1);
			a1[0]=dat1(p1a,0);//obj->getcoordinate(p1a,0);
			a1[1]=dat1(p1a,1);//obj->getcoordinate(p1a,1);
			a2[0]=dat1(p1b,0);//obj->getcoordinate(p1b,0);
			a2[1]=dat1(p1b,1);//obj->getcoordinate(p1b,1);


			vector1<double> uv3(dimension),uv4(dimension);


			double dis3,dis4;

			(obj->getgeo()).distance_vector(a2,a1,uv3,dis3);

			vector1<double> tempvector = a1+(t*(uv3));
			(obj->getgeo()).correct_position(tempvector);
			(obj->getgeo()).distance_vector(a0,tempvector,uv4,dis4);
			//double tt3;
			//this->minimal_distance_between_point_and_line(a0, a1, a2, dis4, uv4, tt3);
			//double t = tt3+0.5;
			double f1 = 0;
			
			f1 = (*bindp).force(sqrt(dis4));
			

			
			//when tt = -0.5 = a1;
			//when tt = 0.5 = a2;
			//matrix<double> forces(3,dimension);

			//tt goes from -0.5 to 0.5

			for(int j = 0 ; j < dimension ; j++) {
			forces(p2,j)+=f1*uv4[j]/sqrt(dis4);
			forces(p1a,j)+=-(1-t)*f1*uv4[j]/sqrt(dis4);
			forces(p1b,j)+=-(t)*f1*uv4[j]/sqrt(dis4);
			}

				cout << "tt: " << tt << endl;
				cout << "t: " <<  t << endl;
				cout << "p2: " << p2 << endl;
				cout << "p1a: " << p1a << endl;
				cout << "p1b: " << p1b << endl;
				cout << "a0: " << a0 << endl;
				cout << "a1: " << a1 << endl;
				cout << "a2: " << a2 << endl;
				cout << "f1: " << f1 << endl;
				cout << "dis4: " << sqrt(dis4) << endl;
				cout << "uv4: " << uv4 << endl;
				cout << endl;


			// if(chckmatrixsize(forces,1000)) {
			// 	cout << "tt: " << tt << endl;
			// 	cout << "t: " <<  t << endl;
			// 	cout << "p2: " << p2 << endl;
			// 	cout << "p1a: " << p1a << endl;
			// 	cout << "p1b: " << p1b << endl;
			// 	cout << "a0: " << a0 << endl;
			// 	cout << "a1: " << a1 << endl;
			// 	cout << "a2: " << a2 << endl;
			// 	cout << "f1: " << f1 << endl;
			// 	cout << "dis4: " << sqrt(dis4) << endl;
			// 	cout << "uv4: " << uv4 << endl;
			// 	pausel();
			// }



		}
	}

	return forces;
}

matrix<double> Microtubule::BindingForces() {
		matrix<double> forces(totalN,dimension);

		for(int i = 0 ; i < na+nb ; i++ ) {
		int p2 = i;
		int p1 = bound[i];
		if(p1>0) {

			double tt  =  bound_along[i];
			int mt =  bound[i];
			double tt2 = tt*(double)(L-1);


			
			int k2 = ceil(tt2);
			int k1 = k2 - 1;

			double t = tt2-k1;

			//na+nb+mt
			int p1a = na+nb+(mt-1)*L+k1;
			int p1b = na+nb+(mt-1)*L+k2;


			vector1<double> a0(dimension),a1(dimension),a2(dimension);

			a0[0]=obj->getcoordinate(p2,0);
			a0[1]=obj->getcoordinate(p2,1);
			a1[0]=obj->getcoordinate(p1a,0);
			a1[1]=obj->getcoordinate(p1a,1);
			a2[0]=obj->getcoordinate(p1b,0);
			a2[1]=obj->getcoordinate(p1b,1);



			vector1<double> uv3(dimension),uv4(dimension);


			double dis3,dis4;

			(obj->getgeo()).distance_vector(a2,a1,uv3,dis3);

			vector1<double> tempvector = a1+(t*(uv3));
			(obj->getgeo()).correct_position(tempvector);

			(obj->getgeo()).distance_vector(a0,tempvector,uv4,dis4);
			//double tt3;
			//this->minimal_distance_between_point_and_line(a0, a1, a2, dis4, uv4, tt3);
			//double t = tt3+0.5;
			double f1 = 0;
			
			f1 = (*bindp).force(sqrt(dis4));
			

			
			//when tt = -0.5 = a1;
			//when tt = 0.5 = a2;
			//matrix<double> forces(3,dimension);

			//tt goes from -0.5 to 0.5

			for(int j = 0 ; j < dimension ; j++) {
			forces(p2,j)+=f1*uv4[j]/sqrt(dis4);

			//BINDING FORCES DON"T ACT ON MICTOTUBULE
			forces(p1a,j)+=-(1-t)*f1*uv4[j]/sqrt(dis4);
			forces(p1b,j)+=-(t)*f1*uv4[j]/sqrt(dis4);
			}
			// cout <<= f1*uv4/sqrt(dis4);
			// cout << "\n";
			// cout <<= -(1-t)*f1*uv4/sqrt(dis4);
			// cout << "\n";
			// cout <<= -(t)*f1*uv4/sqrt(dis4);
			// cout << "\n";
			// if(true) {
			// 	cout << "tt: " << tt << endl;
			// 	cout << "t: " <<  t << endl;
			// 	cout << "p2: " << p2 << endl;
			// 	cout << "p1a: " << p1a << endl;
			// 	cout << "p1b: " << p1b << endl;
			// 	cout << "a0: " << a0 << endl;
			// 	cout << "a1: " << a1 << endl;
			// 	cout << "a2: " << a2 << endl;
			// 	cout << "f1: " << f1 << endl;
			// 	cout << "dis4: " << sqrt(dis4) << endl;
			// 	cout << "uv4: " << uv4 << endl;
			// 	cout << "fp2:" << f1*uv4/sqrt(dis4) << endl;
			// 	cout << "fp1a:" << -(1-t)*f1*uv4/sqrt(dis4) << endl;
			// 	cout << "fp1ab" << -(t)*f1*uv4/sqrt(dis4) << endl;
			// 	pausel();
			// }



		}
	}
	// vector1<double> tbf(2);
	// for(int i = na+nb ; i < na+nb+nc ; i++) {
	// tbf[0] += forces(i,0);
	// tbf[1] += forces(i,1);
	// }
	// // vector1<double> tbfp(2);
	// // for(int i = 0 ; i < na+nb ; i++) {
	// // tbfp[0] += forces(i,0);
	// // tbfp[1] += forces(i,1);
	// // }	
	// cout << tbf[0] << "," << tbf[1] << endl;
	// cout << tbfp[0] << "," << tbfp[1] << endl;
	// pausel();

	return forces;
}

/*
matrix<double> Microtubule::BindingForces() {
	matrix<double> posforces(totalN,dimension);

	for(int i = 0 ; i < na+nb ; i++ ) {
		int p2 = i;
		int p1 = bound[i];

		if(p1>0) {
		double theta1 = (*orient)(p1-na-nb,0);


		vector1<double> a1(dimension);
		vector1<double> a2(dimension);

		vector1<double> a0(dimension);

		a0[0]=obj->getcoordinate(p2,0);
		a0[1]=obj->getcoordinate(p2,1);

		a1[0]=obj->getcoordinate(p1,0)+(L/2.)*cos(theta1);
		a1[1]=obj->getcoordinate(p1,1)+(L/2.)*sin(theta1);

		a2[0]=obj->getcoordinate(p1,0)-(L/2.)*cos(theta1);
		a2[1]=obj->getcoordinate(p1,1)-(L/2.)*sin(theta1);	

		(obj->getgeo()).correct_position(a1);
		(obj->getgeo()).correct_position(a2);

		vector1<double> uv4(dimension);
		double dis4 = 0.0;
		double tt;





		

		//centre of mass of the tube
		vector1<double> c1(dimension);
		c1[0]=obj->getcoordinate(p1,0);
		c1[1]=obj->getcoordinate(p1,1);


 	    this->minimal_distance_between_point_and_line(a0, a1, a2, dis4, uv4, tt);

		double f1 = (*bindp).force(sqrt(dis4));
		for(int j = 0 ; j < dimension ; j++) {
			double fac2 = f1*uv4[j]/sqrt(dis4);
			//posforces(p1,j) += fac2; //BINDING FORCE DOES NOT MOVE MICROTUBULE
			posforces(p2,j) += -fac2; //binding force moves
			}

		}



	}
	return posforces;

}
*/


matrix<double> Microtubule::TubeForces(int &p2,int &p1a, int &p1b) {
	// int p1a = particles[start].b;
	// int p1b = particles[start+1].b;
	vector1<double> a0(dimension),a1(dimension),a2(dimension);
	a0[0]=obj->getcoordinate(p2,0);
	a0[1]=obj->getcoordinate(p2,1);
	a1[0]=obj->getcoordinate(p1a,0);
	a1[1]=obj->getcoordinate(p1a,1);
	a2[0]=obj->getcoordinate(p1b,0);
	a2[1]=obj->getcoordinate(p1b,1);


	vector1<double> uv4(dimension);
	double dis4;
	double tt;
	this->minimal_distance_between_point_and_line(a0, a1, a2, dis4, uv4, tt);
	double t = tt+0.5;
	double f1 = 0;
	if(p2<na) {
		f1 = (*fac).force(sqrt(dis4));
	}
	else {
		f1 = (*fbc).force(sqrt(dis4));
	}
	//when tt = -0.5 = a1;
	//when tt = 0.5 = a2;
	matrix<double> forces(3,dimension);

	//tt goes from -0.5 to 0.5

	for(int j = 0 ; j < dimension ; j++) {
	forces(0,j)=-f1*uv4[j]/sqrt(dis4);
	forces(1,j)=(1-t)*f1*uv4[j]/sqrt(dis4);
	forces(2,j)=(t)*f1*uv4[j]/sqrt(dis4);
	}

	return forces;

}

void Microtubule::ForcesDueToPositionPL(matrix<int> &pairs, matrix<double> &forces) {
	//we have pairs, which is a list of free particles 
	if(pairs.getNsafe() == 0) return;
	vector<mdpair> particles;
	for(int i = 0 ; i < pairs.getNsafe() ; i++ ) {
		int p1 = pairs(i,0) > pairs(i,1) ?  pairs(i,0) : pairs(i,1) ;
		int p2 = pairs(i,0) < pairs(i,1) ?  pairs(i,0) : pairs(i,1) ;
		mdpair pl(p2,p1);
		particles.push_back(pl);
	}
	sort(particles.begin(),particles.end()); //SORT THE LIST WITH RESPECT TO THE FREE PARTICLES
	//cout << "sorted " << pairs.getNsafe() << endl;
	// cout << particles.size() << endl;
	// for(int i = 0 ; i < particles.size() ; i++ ) {
	// 	cout << particles[i].a <<  " " << particles[i].b << endl;
	// }
	// cout << endl;
	// pausel();
	//from this list, find the continual runs of
	int iter  = 0;
	int start = iter;
	int kl1 = particles[iter].a;
	int mt1 = (int)((double)((particles[iter]).b-na-nb)/(double)(L));
	for(;;) {
		//find the point which corresponds to the interaction of the particle with a single microtubule
		iter++;
		// cout << iter << endl;
		// pausel();
		//cout << "p2: " << iter << endl;
		if(iter>particles.size() ) break;
		int kl2 = particles[iter].a;
		int mt2 = (int)((double)((particles[iter]).b-na-nb)/(double)(L));

		if(kl1!=kl2 || mt1 != mt2) {
			int finish = iter;
			mt1 = mt2;
			kl1 = kl2;

			if(finish-start == 0) {
				error("something wrong");
			}
			else if(finish-start == 1 ) { //if it is close to only one of the points
				int p2 = particles[start].a;
				int p1 = particles[start].b;


				vector1<double> uv4(dimension);
				double dis4;
				obj->disvec(p1,p2,uv4,dis4);

				double f1 = 0;
				if(p2<na) {
					f1 = (*fac).force(sqrt(dis4));
				}
				else {
					f1 = (*fbc).force(sqrt(dis4));
				}



				for(int j = 0 ; j < dimension ; j++) {
					double fac2 = f1*uv4[j]/sqrt(dis4);

					forces(p1,j) += fac2;
					forces(p2,j) += -fac2;

				}
				//cout << start <<  " " << finish << endl;

			}
			else if(finish-start == 2) { //if it is close to two of the points


				int p2 = particles[start].a;

				int p1a = particles[start].b;
				int p1b = particles[start+1].b;

				if(p1b < p1a) {
					int p1at = p1a;
					int p1bt = p1b;
					p1a =  p1bt;
					p1b =  p1at;

				}

				matrix<double> tempforce = this->TubeForces(p2,p1a,p1b);


				for(int j = 0 ; j < dimension ; j++) {
				forces(p2,j)+=tempforce(0,j);
				forces(p1a,j)+=tempforce(1,j);
				forces(p1b,j)+=tempforce(2,j);
				}
			//	cout << start <<  " " << finish << endl;


			}
			else { //if it close to more of the points
				vector<dispair> distances;
				int p2 = particles[start].a;
				for(int i = start ; i < finish ; i++ ) {
					int p1 = particles[i].b;
					//int p2 = particles[i].a;

					vector1<double> uv(dimension);
					double dissqr;
					obj->disvec(p1,p2,uv,dissqr);
					dispair hg(p1,dissqr);

					distances.push_back(hg);
				}
				sort(distances.begin(),distances.end());

				int p1a = distances[0].a;
				int p1b = distances[1].a;

				if(p1b < p1a) {
					int p1at = p1a;
					int p1bt = p1b;
					p1a =  p1bt;
					p1b =  p1at;

				}

				matrix<double> tempforce = this->TubeForces(p2,p1a,p1b);


				for(int j = 0 ; j < dimension ; j++) {
				forces(p2,j)+=tempforce(0,j);
				forces(p1a,j)+=tempforce(1,j);
				forces(p1b,j)+=tempforce(2,j);
				}
			//	cout << start <<  " " << finish << endl;


			}

			start = finish;
		}

	}

//	pausel();

	//sorted list
}



// void Microtubule::UpdateBoundAlong() {
//  	for(int i = 0  ; i < na + nb ; i++ ) {
//  		if(bound[i]>0) {


// 				vector<dispair> distances;
// 				int p2 = i;
// 				for(int j = 0 ; j < L ; j++ ) {
// 					int p1 =na+nb+(mt-1)*L+j
// 					//int p1 = particles[i].b;
// 					//int p2 = particles[i].a;

// 					vector1<double> uv(dimension);
// 					double dissqr;
// 					obj->disvec(p1,p2,uv,dissqr);
// 					dispair hg(p1,dissqr);

// 					distances.push_back(hg);
// 				}
// 				sort(distances.begin(),distances.end());

// 				int p1a = distances[0].a;
// 				int p1b = distances[1].a;
// 				if( p1b < p1a ) {

// 				}






//  		}
//  	}
// }

