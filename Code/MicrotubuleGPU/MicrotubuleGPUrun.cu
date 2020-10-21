
void Microtubule::runGPU(int runtime, int every)
{
//	pausel();
	//int ccc;
	int totalN = obj->getN();

	//num is the number of boxes per length

	int ncells = num*num;

	WCApotentialGPU faa_gpu(2.,1.,2.);
	WCApotentialGPU fab_gpu(1.,1.,0.);
	WCApotentialGPU fac_gpu(1.,1.,0.);
	WCApotentialGPU fbb_gpu(2.,1.,2.);
	WCApotentialGPU fbc_gpu(1.,1.,0.);
	WCApotentialGPU fcc_gpu(1.,1.,0.);
	HarmonicPotentialGPU bindp_gpu(100.,1.); 
	FENEPotentialGPU bindm_gpu(50.,1.5); 
	BendingPotentialGPU bendp_gpu(100.,0.);

//we now have the count of each cell list
	int nbpairs = 5*ncells;
	int nperl = num;

	int *cells1 = new int [nbpairs];
	int *cells2 = new int [nbpairs];


	int itery = 0;
	for(int i1 = 0 ; i1 < num ; i1++) {
		for(int i2 = 0 ; i2 < num ; i2++ ) {


			int b1 =  i1*nperl+i2;

			int i3 = i1+0;
			int j3 = i2+0;

			int i4 = i1+1;
			int j4 = i2+0;

			int i5 = i1-1;
			int j5 = i2+1;

			int i6 = i1+0;
			int j6 = i2+1;

			int i7 = i1+1;
			int j7 = i2+1;

			prdshft(i3,nperl);
			prdshft(j3,nperl);

			prdshft(i4,nperl);
			prdshft(j4,nperl);

			prdshft(i5,nperl);
			prdshft(j5,nperl);
			
			prdshft(i6,nperl);
			prdshft(j6,nperl);
			
			prdshft(i7,nperl);
			prdshft(j7,nperl);		

			cells1[itery] =  b1;
			cells2[itery] =  i3*nperl+j3;

			itery++;

			cells1[itery] =  b1;
			cells2[itery] =  i4*nperl+j4;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i5*nperl+j5;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i6*nperl+j6;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i7*nperl+j7;

			itery++;


		}
	}
	int size4 = nbpairs*sizeof(int);

	int *d_cells1;
	int *d_cells2;

	cudaMalloc((void**)&d_cells1,size4);

	cudaMalloc((void**)&d_cells2,size4);

	cudaMemcpy(d_cells1,cells1,size4,cudaMemcpyHostToDevice);
	cudaMemcpy(d_cells2,cells2,size4,cudaMemcpyHostToDevice);
	//int sibdiv = floor(ll/4.0);
	// print_device_array(d_cells1,nbpairs);
	// print_device_array(d_cells2,nbpairs);



	// matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);
	
	float2 *particles = new float2 [totalN];
	float2 *momenta = new float2 [totalN];
	int *p_indices = new int [totalN];

	for(int i = 0 ; i < totalN ; i++)
	p_indices[i]=i;

	float2 *d_particles;
	float2 *d_momenta;
	int *d_p_indices;

	int *d_bound;
	double *d_boundalong;
	int *d_changestate;

	cudaMalloc((void**)&d_bound,(na+nb)*sizeof(int));
	cudaMalloc((void**)&d_boundalong,(na+nb)*sizeof(double));
	cudaMalloc((void**)&d_changestate,(na+nb)*sizeof(int));

	cudaMemset(d_bound,0,(na+nb)*sizeof(int));
	cudaMemset(d_boundalong,0.,(na+nb)*sizeof(double));
	cudaMemset(d_changestate,0,(na+nb)*sizeof(int));

	double *d_totalforcex;
	double *d_totalforcey;

	cudaMalloc((void**)&d_totalforcex,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey,totalN*sizeof(double));


	matrix<double> state(obj->getdat());


	for(int i = 0  ; i < totalN ; i++) {

	float2 c;
	c.x=state(i,0);
	c.y=state(i,1);

	(particles)[i]=c;

	float2 d;

	d.x = 0.;
	d.y = 0.;

	(momenta)[i]=d;
	}


	int size =  totalN*sizeof(float2);
	int size2 = totalN*sizeof(int);


	cudaMalloc((void**)&d_particles,size);
	cudaMalloc((void**)&d_momenta,size);
	cudaMalloc((void**)&d_p_indices,size2);

	cudaMemcpy(d_particles,particles,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_momenta,momenta,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_p_indices,p_indices,size2,cudaMemcpyHostToDevice);


	// matrix<int> *froyo1 = obj->calculatepairs(boxes,pai,3.5);
	// matrix<int> *froyo2 = obj->calculatepairs(boxes,pbi,3.5);
	// matrix<int> *froyo3 = obj->calculatepairs(boxes,pci,3.5);
	// matrix<int> *froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
	// matrix<int> *froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
	// matrix<int> *froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);
	
	int *d_indices1;
	int *d_indices2;
	double *d_close;


	int tpp;
	


	construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp);



	less_than_condition_AND cond1(SQR(3.5),0,na);
	less_than_condition_AND cond2(SQR(3.5),na,na+nb);
	less_than_condition_AND cond3(SQR(3.5),na+nb,na+nb+nc);
	less_than_condition_NAND cond4(SQR(3.5),0,na,na,na+nb);
	less_than_condition_NAND cond5(SQR(3.5),0,na,na+nb,na+nb+nc);
	less_than_condition_NAND cond6(SQR(3.5),na,na+nb,na+nb,na+nb+nc);


	int th1;
	int *d1_list1,*d1_list2,*d1_list3,*d1_list4;
	pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa

	int th2;
	int *d2_list1,*d2_list2,*d2_list3,*d2_list4;
	pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2);	//fbb

	int th3;
	int *d3_list1,*d3_list2,*d3_list3,*d3_list4;
	pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3);	//fcc

	int th4;
	int *d4_list1,*d4_list2,*d4_list3,*d4_list4;
	pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4);	//fab

	int th5;
	int *d5_list1,*d5_list2,*d5_list3,*d5_list4;
	pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5);	//fac

	int th6;
	int *d6_list1,*d6_list2,*d6_list3,*d6_list4;
	pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc			
	//matrix<double> state(obj->getdat()); //the state of the system

	int th8 = (*bondpairs).getNsafe();
	int th9 = (*bendtriplets).getNsafe();


	int *d8_list1,*d8_list2,*d8_list3,*d8_list4;
	cudaMalloc((void**)&d8_list1,th8*sizeof(int));
	cudaMalloc((void**)&d8_list2,th8*sizeof(int));
	cudaMalloc((void**)&d8_list3,th8*sizeof(int));
	cudaMalloc((void**)&d8_list4,th8*sizeof(int));



	int *h8_list1 = new int [th8];
	int *h8_list2 = new int [th8];

	for(int i = 0 ; i < th8 ; i++ ) {
		h8_list1[i] = (*bondpairs)(i,0);
		h8_list2[i] = (*bondpairs)(i,1);
	}

	cudaMemcpy(d8_list1,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list2,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list3,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list4,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);





	int *d9_list1,*d9_list2,*d9_list3,*d9_list4,*d9_list5,*d9_list6;
	cudaMalloc((void**)&d9_list1,th9*sizeof(int));
	cudaMalloc((void**)&d9_list2,th9*sizeof(int));
	cudaMalloc((void**)&d9_list3,th9*sizeof(int));
	cudaMalloc((void**)&d9_list4,th9*sizeof(int));
	cudaMalloc((void**)&d9_list5,th9*sizeof(int));
	cudaMalloc((void**)&d9_list6,th9*sizeof(int));

	int *h9_list1 = new int [th9];
	int *h9_list2 = new int [th9];
	int *h9_list3 = new int [th9];

	for(int i = 0 ; i < th9 ; i++ ) {
		h9_list1[i] = (*bendtriplets)(i,0);
		h9_list2[i] = (*bendtriplets)(i,1);
		h9_list3[i] = (*bendtriplets)(i,2);
	}	

	cudaMemcpy(d9_list1,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list2,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list3,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list4,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list5,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list6,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);


	int i;

	double cons1;
	double cons2;
	double cons3;
	double cons4;

	//(mom)->operator()(i,i1) = c5*c2*((mom)->operator()(i,i1)) + (c5*(c3)+q)*F(i,i1) + (c5*(c4)+r)*R(i,i1);
	cons1 = (*obj).getc5()*(*obj).getc2();
	cons2 = (*obj).getc5()*(*obj).getc3()+(*obj).getq();
	cons3 = (*obj).getc5()*(*obj).getc4()+(*obj).getr();
	cons4 = (*obj).getc1();

	double d_dt = (*obj).getdt();
	double d_m = (*obj).getm();
	double d_kT = (*obj).getkT();


	// vector<matrix<double> > savef1forces;
	// vector<matrix<double> > savef2forces;
	// vector<matrix<double> > savef3forces;
	// vector<matrix<double> > savef4forces;
	// vector<matrix<double> > savef5forces;
	// vector<matrix<double> > savef6forces;
	// vector<matrix<double> > savef7forces;

	// vector<matrix<double> > saveftemp1forces;
	// vector<matrix<double> > saveftemp2forces;
	// vector<matrix<double> > saveftemp3forces;

	// vector<matrix<double> > savepositions;
	// vector<vector1<int> > savebound;
	// vector<vector1<double> > saveboundalongs;

	for(i = 0 ; i < runtime ; i++) {
		//cout << i << endl;

		cout << i << endl;
	
		//cout << (*obj).avmom() << endl;
	if(i%25==0) {
		//delete froyo1,froyo2,froyo3,froyo4,froyo5,froyo6;
		// cout << "updated after: " << i << endl;
		// state = obj->getdat();


		cudaFree(d1_list1);cudaFree(d1_list2);cudaFree(d1_list3);cudaFree(d1_list4);
		

		cudaFree(d2_list1);cudaFree(d2_list2);cudaFree(d2_list3);cudaFree(d2_list4);
		

		cudaFree(d3_list1);cudaFree(d3_list2);cudaFree(d3_list3);cudaFree(d3_list4);
		

		cudaFree(d4_list1);cudaFree(d4_list2);cudaFree(d4_list3);cudaFree(d4_list4);
		

		cudaFree(d5_list1);cudaFree(d5_list2);cudaFree(d5_list3);cudaFree(d5_list4);
	
		cudaFree(d6_list1);cudaFree(d6_list2);cudaFree(d6_list3);cudaFree(d6_list4);

		cudaFree(d_indices1);

		cudaFree(d_indices2);

		cudaFree(d_close);

		this->resetindices(d_p_indices,totalN);


		construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp,false);

		// cout << tpp << endl;
		// cout << "pair" << endl;

		pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa
		pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2); //fbb
		pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3); //fcc
		pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4); //fab
		pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5); //fac
		pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc

	}




	double *d1_forces1x;
	double *d1_forces2x;
	double *d1_forces1y;
	double *d1_forces2y;
	cudaMalloc((void**)&d1_forces1x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces1y,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2y,th1*sizeof(double));
	calculateforces2D(d1_list1,d1_list2,d_particles, d1_forces1x,d1_forces1y,d1_forces2x,d1_forces2y, faa_gpu ,th1, l,true);

	// arracychck(d1_forces1x,th1);
	// arracychck(d1_forces1y,th1);
	// cout << "force1" << endl;

	double *d2_forces1x;
	double *d2_forces2x;
	double *d2_forces1y;
	double *d2_forces2y;
	cudaMalloc((void**)&d2_forces1x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces1y,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2y,th2*sizeof(double));
	calculateforces2D(d2_list1,d2_list2,d_particles, d2_forces1x,d2_forces1y,d2_forces2x,d2_forces2y, fbb_gpu ,th2, l,true);	



	// arracychck(d2_forces1x,th2);
	// arracychck(d2_forces1y,th2);

	// cout << "force2" << endl;
	double *d3_forces1x;
	double *d3_forces2x;
	double *d3_forces1y;
	double *d3_forces2y;
	cudaMalloc((void**)&d3_forces1x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces1y,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2y,th3*sizeof(double));
	calculateforces2D(d3_list1,d3_list2,d_particles, d3_forces1x,d3_forces1y,d3_forces2x,d3_forces2y, fcc_gpu ,th3, l,true);	

	// arracychck(d3_forces1x,th3);
	// arracychck(d3_forces1y,th3);

	// cout << "force3" << endl;
	double *d4_forces1x;
	double *d4_forces2x;
	double *d4_forces1y;
	double *d4_forces2y;
	cudaMalloc((void**)&d4_forces1x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces1y,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2y,th4*sizeof(double));
	calculateforces2D(d4_list1,d4_list2,d_particles, d4_forces1x,d4_forces1y,d4_forces2x,d4_forces2y, fab_gpu ,th4, l,true);
	
	// arracychck(d4_forces1x,th4);
	// arracychck(d4_forces1y,th4);
	// cout << "force4" << endl;

	double *d5_forces1x;
	double *d5_forces2x;
	double *d5_forces1y;
	double *d5_forces2y;
	cudaMalloc((void**)&d5_forces1x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces1y,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2y,th5*sizeof(double));
	calculateforces2D(d5_list1,d5_list2,d_particles, d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y, fbc_gpu ,th5, l,true);
	// arracychck(d5_forces1x,th5);
	// arracychck(d5_forces1y,th5);
	// cout << "force5" << endl;
	// cout << "d5" << endl;
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y,d_particles,th5,totalN);
	// pausel();	

	double *d6_forces1x;
	double *d6_forces2x;
	double *d6_forces1y;
	double *d6_forces2y;
	cudaMalloc((void**)&d6_forces1x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces1y,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2y,th6*sizeof(double));
	calculateforces2D(d6_list1,d6_list2,d_particles, d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y, fac_gpu ,th6, l,true);	
	// arracychck(d6_forces1x,th6);
	// arracychck(d6_forces1y,th6);
	// cout << "force6" << endl;
	// cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();

	// matrix<double> F6((*obj).calculateforces(*bondpairs,*bindm));
	callCalculateUnbindingsGPU(d_particles,d_bound,d_boundalong,d_changestate);



	// cout << "unbindings calculated" << endl;


	callCalculateBindingsGPU(d5_list1,d5_list2,d6_list1,d6_list2,d_particles, d_bound, d_boundalong,d_changestate,th5 ,th6 );


	// cout << "bindings calculated" << endl;


	int *d7_list1,*d7_list2,*d7_list3;
	double *d7_forces1x;
	double *d7_forces1y;
	double *d7_forces2x;
	double *d7_forces2y;
	double *d7_forces3x;
	double *d7_forces3y;
	int th7;
	BindingForcesGPU(d_particles, d_bound, d_boundalong, d7_list1,d7_list2,d7_list3, d7_forces1x, d7_forces1y, d7_forces2x,d7_forces2y,d7_forces3x, d7_forces3y, bindp_gpu, th7);

	// arracychck(d7_forces1x,th7);
	// arracychck(d7_forces1y,th7);
	// cout << "binding forces calculated" << endl;


	double *d8_forces1x;
	double *d8_forces2x;
	double *d8_forces1y;
	double *d8_forces2y;
	cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));
	calculateforces2D(d8_list1,d8_list2,d_particles, d8_forces1x,d8_forces1y,d8_forces2x,d8_forces2y, bindm_gpu ,th8, l,true);
	
	// arracychck(d8_forces1x,th8);
	// arracychck(d8_forces1y,th8);

	// cout << "force 8" << endl;
	resetchangestate(d_changestate);
	// matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets,*bendp));

//	int th9;


	double *d9_forces1x;
	double *d9_forces2x;
	double *d9_forces1y;
	double *d9_forces2y;
	double *d9_forces3x;
	double *d9_forces3y;
	// cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));	
	cudaMalloc((void**)&d9_forces1x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces1y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3y,th9*sizeof(double));
	BendingForcesGPU(d_particles, d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces1y,d9_forces2x,d9_forces2y, d9_forces3x, d9_forces3y,bendp_gpu,th9);

	// print_device_array_weave(d9_forces1x,d9_forces1y,th9);
	// print_device_array_weave(d9_forces2x,d9_forces2y,th9);
	// print_device_array_weave(d9_forces3x,d9_forces3y,th9);
	// cout << "force 9" << endl;	
	// arracychck(d9_forces1x,th9);
	// arracychck(d9_forces1y,th9);

	// print_device_float2(d_particles,totalN);
	// print_device_array(d9_list1,th9);
	// print_device_array(d9_list2,th9);
	// print_device_array(d9_list3,th9);

	// cout << "force9" << endl;

	//	matrix<double> F = ftemp1+ftemp2+ftemp3+F1+F2+F3+F4+F5+F6+F7;//+F4+F5;

		// matrix<double> R(totalN,dimension);
		// for(int i1 = 0 ; i1 < totalN ; i1++) {
		// 	for(int j = 0 ; j < dimension ; j++) {
		// 		R(i1,j) = (3.464101615 * ((double) rand() / (RAND_MAX)) - 1.732050808);
		// 	}
		// }

	int *d10_list1;
	double *d10_forces1x;
	double *d10_forces1y;	
	int th10;

	PositionForcesDueToAnglesGPU(d_particles, d_bound, d_boundalong, d10_list1, d10_forces1x, d10_forces1y,th10);
	

	// arracychck(d10_forces1x,th10);
	// arracychck(d10_forces1y,th10);
	// cout << "force10" << endl;
	// cout << "all forces calculated" << endl;





	resetforce(d_totalforcex);

	resetforce(d_totalforcey);

	//cout << "reset" << endl;

	// print_device_array(d_totalforcex,totalN);


	ReduceForces(d1_list1,d1_list2,d1_list3,d1_list4,d1_forces1x,d1_forces2x,d1_forces1y,d1_forces2y,d_totalforcex,d_totalforcey,th1);
	// cout << "d1" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d2_list1,d2_list2,d2_list3,d2_list4,d2_forces1x,d2_forces2x,d2_forces1y,d2_forces2y,d_totalforcex,d_totalforcey,th2);
	/// cout << "d2" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces(d3_list1,d3_list2,d3_list3,d3_list4,d3_forces1x,d3_forces2x,d3_forces1y,d3_forces2y,d_totalforcex,d_totalforcey,th3);
	// cout << "d3" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d4_list1,d4_list2,d4_list3,d4_list4,d4_forces1x,d4_forces2x,d4_forces1y,d4_forces2y,d_totalforcex,d_totalforcey,th4);
	// cout << "d4" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d5_list1,d5_list2,d5_list3,d5_list4,d5_forces1x,d5_forces2x,d5_forces1y,d5_forces2y,d_totalforcex,d_totalforcey,th5);
	// cout << "d5" << endl;	
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d_particles,th5,totalN);
	// pausel();	
	ReduceForces(d6_list1,d6_list2,d6_list3,d6_list4,d6_forces1x,d6_forces2x,d6_forces1y,d6_forces2y,d_totalforcex,d_totalforcey,th6);	
	// cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();
	ReduceForces3(d7_list1,d7_list2,d7_list3,d7_forces1x,d7_forces2x,d7_forces3x,d7_forces1y,d7_forces2y,d7_forces3y,d_totalforcex,d_totalforcey,th7);	
	// cout << "d7" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces(d8_list1,d8_list2,d8_forces1x,d8_forces2x,d8_forces1y,d8_forces2y,d_totalforcex,d_totalforcey,th8);
	// cout << "d8" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces3(d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces2x,d9_forces3x,d9_forces1y,d9_forces2y,d9_forces3y,d_totalforcex,d_totalforcey,th9);		
	
	// cout << "d9" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double ff = (v0_a+v0_b)/2.;
	ReduceForcesAndNormalize(d10_list1,d10_forces1x,d10_forces1y,d_totalforcex,d_totalforcey, max_s, ff, th10);

	cout << "reduction" << endl;
	// cout << "d10" << endl;
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	if(i>0&&i%every==0) { 
		// for(int j = 0 ; j < na+nb ; j++) {
		// if(bound[j]>0){
		// cout << "printed" << endl;
		// cout << j << endl;
		// cout << bound[j] << endl;
		// cout << bound_along[j] << endl;
		// cout << F5(j,'r') << endl;
		// cout << F4(j,'r') << endl;
		// cout << obj->getcoordinate(j,0) << " " << obj->getcoordinate(j,1) << endl;
		// cout << obj->getcoordinate(200,0) <<  " " << obj->getcoordinate(200,1) << endl;
		// cout << obj->getcoordinate(200+L,0) <<  " " << obj->getcoordinate(200+L,1) << endl;
		// }
		// }
		// cout << F6 << endl;
		// cout << F7 << endl;

		// cout << ftemp2 << endl;
		// cout << ftemp3 << endl;

		stringstream ss2;
		// ss2 <<i/every;
		// string pairlist = "list";

		 stringstream kts;
		 kts << (*obj).getkT();

		 // stringstream epi;
		 // epi << eps;

		 // stringstream epieq;
		 // epieq << eqeps;

		 stringstream len;
		 len << l;

		 string extension =  "_kT="+kts.str()+"_l="+len.str()+".csv";

		stringstream ss;
		ss <<(i/every);
		string filename = "x";
		filename += ss.str();
		filename += extension;

		string momname = "bind";
		momname += ss.str();
		momname += extension;

		string baname = "bind_along";
		baname += ss.str();
		baname += extension;



		ofstream myfile;
		myfile.open(filename.c_str());
		//myfile <<= (*obj).getdat();
		file_print_device_float2(d_particles,totalN,myfile);
		myfile.close();



		}	

	double *d_R1;
	double *d_R2;

	// cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));	
	cudaMalloc((void**)&d_R1,totalN*sizeof(double));
	cudaMalloc((void**)&d_R2,totalN*sizeof(double));

	setstaterandom(d_R1,1.732050808,totalN);
 	setstaterandom(d_R2,1.732050808,totalN);	

 	//advmom2D(d_momenta, d_totalforcex, d_totalforcey, d_R1, d_R2, cons1,cons2,cons3,totalN);
// 	advmom2_spatialdependence(d_momenta,d_particles,d_totalforcex,d_totalforcey,d_R1,d_R2,func, d_dt, d_kT, d_m, totalN);
	advmom2D(d_momenta,d_totalforcex,d_totalforcey,d_R1,d_R2,cons1,cons2,cons3, totalN);
 	advpos2D(d_particles, d_momenta, cons4, totalN);


 	applypbc2D(d_particles,d_momenta,l,is_periodic,totalN);

 	cout << "updated" << endl;

	cudaFree(d1_forces1x);
	cudaFree(d1_forces2x);
	cudaFree(d1_forces1y);
	cudaFree(d1_forces2y);
	cudaFree(d2_forces1x);
	cudaFree(d2_forces2x);
	cudaFree(d2_forces1y);
	cudaFree(d2_forces2y);
	cudaFree(d3_forces1x);
	cudaFree(d3_forces2x);
	cudaFree(d3_forces1y);
	cudaFree(d3_forces2y);
	cudaFree(d4_forces1x);
	cudaFree(d4_forces2x);
	cudaFree(d4_forces1y);
	cudaFree(d4_forces2y);
	cudaFree(d5_forces1x);
	cudaFree(d5_forces2x);
	cudaFree(d5_forces1y);
	cudaFree(d5_forces2y);
	cudaFree(d6_forces1x);
	cudaFree(d6_forces2x);
	cudaFree(d6_forces1y);
	cudaFree(d6_forces2y);
 	cudaFree(d7_list1);
 	cudaFree(d7_list2);
 	cudaFree(d7_list3);
	cudaFree(d7_forces1x);
	cudaFree(d7_forces2x);
	cudaFree(d7_forces3x);
	cudaFree(d7_forces1y);
	cudaFree(d7_forces2y);
	cudaFree(d7_forces3y);
	cudaFree(d8_forces1x);
	cudaFree(d8_forces2x);
	cudaFree(d8_forces1y);
	cudaFree(d8_forces2y);
	cudaFree(d9_forces1x);
	cudaFree(d9_forces2x);
	cudaFree(d9_forces3x);
	cudaFree(d9_forces1y);
	cudaFree(d9_forces2y);
	cudaFree(d9_forces3y);
	cudaFree(d10_list1);
	cudaFree(d10_forces1x);
	cudaFree(d10_forces1y);

	cout << "freed" << endl;


	}
}


/*
template <typename Fun>
void Microtubule::runGPUcheck(int runtime, int every, Fun func)
{
//	pausel();
	int ccc;
	int totalN = obj->getN();

	//num is the number of boxes per length

	int ncells = num*num;

	WCApotentialGPU faa_gpu(2.,1.,2.);
	WCApotentialGPU fab_gpu(1.,1.,0.);
	WCApotentialGPU fac_gpu(1.,1.,0.);
	WCApotentialGPU fbb_gpu(2.,1.,2.);
	WCApotentialGPU fbc_gpu(1.,1.,0.);
	WCApotentialGPU fcc_gpu(1.,1.,0.);
	HarmonicPotentialGPU bindp_gpu(100.,0.); 
	FENEPotentialGPU bindm_gpu(50.,1.5); 
	BendingPotentialGPU bendp_gpu(100.,0.);

//we now have the count of each cell list

	matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);


	matrix<int> *froyo1 = obj->calculatepairs(boxes,pai,3.5);
	matrix<int> *froyo2 = obj->calculatepairs(boxes,pbi,3.5);
	matrix<int> *froyo3 = obj->calculatepairs(boxes,pci,3.5);
	matrix<int> *froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
	matrix<int> *froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
	matrix<int> *froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);


	int nbpairs = 5*ncells;
	int nperl = num;

	int *cells1 = new int [nbpairs];
	int *cells2 = new int [nbpairs];


	int itery = 0;
	for(int i1 = 0 ; i1 < num ; i1++) {
		for(int i2 = 0 ; i2 < num ; i2++ ) {


			int b1 =  i1*nperl+i2;

			int i3 = i1+0;
			int j3 = i2+0;

			int i4 = i1+1;
			int j4 = i2+0;

			int i5 = i1-1;
			int j5 = i2+1;

			int i6 = i1+0;
			int j6 = i2+1;

			int i7 = i1+1;
			int j7 = i2+1;

			prdshft(i3,nperl);
			prdshft(j3,nperl);

			prdshft(i4,nperl);
			prdshft(j4,nperl);

			prdshft(i5,nperl);
			prdshft(j5,nperl);
			
			prdshft(i6,nperl);
			prdshft(j6,nperl);
			
			prdshft(i7,nperl);
			prdshft(j7,nperl);		

			cells1[itery] =  b1;
			cells2[itery] =  i3*nperl+j3;

			itery++;

			cells1[itery] =  b1;
			cells2[itery] =  i4*nperl+j4;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i5*nperl+j5;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i6*nperl+j6;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i7*nperl+j7;

			itery++;


		}
	}
	int size4 = nbpairs*sizeof(int);

	int *d_cells1;
	int *d_cells2;

	cudaMalloc((void**)&d_cells1,size4);

	cudaMalloc((void**)&d_cells2,size4);

	cudaMemcpy(d_cells1,cells1,size4,cudaMemcpyHostToDevice);
	cudaMemcpy(d_cells2,cells2,size4,cudaMemcpyHostToDevice);
	//int sibdiv = floor(ll/4.0);
	// print_device_array(d_cells1,nbpairs);
	// print_device_array(d_cells2,nbpairs);



	// matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);
	
	float2 *particles = new float2 [totalN];
	float2 *momenta = new float2 [totalN];
	int *p_indices = new int [totalN];

	for(int i = 0 ; i < totalN ; i++)
	p_indices[i]=i;

	float2 *d_particles;
	float2 *d_momenta;
	int *d_p_indices;

	int *d_bound;
	double *d_boundalong;
	int *d_changestate;

	cudaMalloc((void**)&d_bound,(na+nb)*sizeof(int));
	cudaMalloc((void**)&d_boundalong,(na+nb)*sizeof(double));
	cudaMalloc((void**)&d_changestate,(na+nb)*sizeof(int));

	cudaMemset(d_bound,0,(na+nb)*sizeof(int));
	cudaMemset(d_boundalong,0.,(na+nb)*sizeof(double));
	cudaMemset(d_changestate,0,(na+nb)*sizeof(int));

	double *d_totalforcex;
	double *d_totalforcey;

	cudaMalloc((void**)&d_totalforcex,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey,totalN*sizeof(double));


	matrix<double> state(obj->getdat());


	for(int i = 0  ; i < totalN ; i++) {

	float2 c;
	c.x=state(i,0);
	c.y=state(i,1);

	(particles)[i]=c;

	float2 d;

	d.x = 0.;
	d.y = 0.;

	(momenta)[i]=d;
	}


	int size =  totalN*sizeof(float2);
	int size2 = totalN*sizeof(int);


	cudaMalloc((void**)&d_particles,size);
	cudaMalloc((void**)&d_momenta,size);
	cudaMalloc((void**)&d_p_indices,size2);

	cudaMemcpy(d_particles,particles,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_momenta,momenta,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_p_indices,p_indices,size2,cudaMemcpyHostToDevice);


	// matrix<int> *froyo1 = obj->calculatepairs(boxes,pai,3.5);
	// matrix<int> *froyo2 = obj->calculatepairs(boxes,pbi,3.5);
	// matrix<int> *froyo3 = obj->calculatepairs(boxes,pci,3.5);
	// matrix<int> *froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
	// matrix<int> *froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
	// matrix<int> *froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);
	
	int *d_indices1;
	int *d_indices2;
	double *d_close;


	int tpp;
	


	construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp);



	less_than_condition_AND cond1(SQR(3.5),0,na);
	less_than_condition_AND cond2(SQR(3.5),na,na+nb);
	less_than_condition_AND cond3(SQR(3.5),na+nb,na+nb+nc);
	less_than_condition_NAND cond4(SQR(3.5),0,na,na,na+nb);
	less_than_condition_NAND cond5(SQR(3.5),0,na,na+nb,na+nb+nc);
	less_than_condition_NAND cond6(SQR(3.5),na,na+nb,na+nb,na+nb+nc);


	int th1;
	int *d1_list1,*d1_list2,*d1_list3,*d1_list4;
	pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa

	int th2;
	int *d2_list1,*d2_list2,*d2_list3,*d2_list4;
	pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2);	//fbb

	int th3;
	int *d3_list1,*d3_list2,*d3_list3,*d3_list4;
	pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3);	//fcc

	int th4;
	int *d4_list1,*d4_list2,*d4_list3,*d4_list4;
	pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4);	//fab

	int th5;
	int *d5_list1,*d5_list2,*d5_list3,*d5_list4;
	pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5);	//fac

	int th6;
	int *d6_list1,*d6_list2,*d6_list3,*d6_list4;
	pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc			
	//matrix<double> state(obj->getdat()); //the state of the system

	int th8 = (*bondpairs).getNsafe();
	int th9 = (*bendtriplets).getNsafe();


	int *d8_list1,*d8_list2,*d8_list3,*d8_list4;
	cudaMalloc((void**)&d8_list1,th8*sizeof(int));
	cudaMalloc((void**)&d8_list2,th8*sizeof(int));
	cudaMalloc((void**)&d8_list3,th8*sizeof(int));
	cudaMalloc((void**)&d8_list4,th8*sizeof(int));



	int *h8_list1 = new int [th8];
	int *h8_list2 = new int [th8];

	for(int i = 0 ; i < th8 ; i++ ) {
		h8_list1[i] = (*bondpairs)(i,0);
		h8_list2[i] = (*bondpairs)(i,1);
	}

	cudaMemcpy(d8_list1,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list2,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list3,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list4,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);





	int *d9_list1,*d9_list2,*d9_list3,*d9_list4,*d9_list5,*d9_list6;
	cudaMalloc((void**)&d9_list1,th9*sizeof(int));
	cudaMalloc((void**)&d9_list2,th9*sizeof(int));
	cudaMalloc((void**)&d9_list3,th9*sizeof(int));
	cudaMalloc((void**)&d9_list4,th9*sizeof(int));
	cudaMalloc((void**)&d9_list5,th9*sizeof(int));
	cudaMalloc((void**)&d9_list6,th9*sizeof(int));

	int *h9_list1 = new int [th9];
	int *h9_list2 = new int [th9];
	int *h9_list3 = new int [th9];

	for(int i = 0 ; i < th9 ; i++ ) {
		h9_list1[i] = (*bendtriplets)(i,0);
		h9_list2[i] = (*bendtriplets)(i,1);
		h9_list3[i] = (*bendtriplets)(i,2);
	}	

	cudaMemcpy(d9_list1,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list2,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list3,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list4,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list5,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list6,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);


	int i;

	double cons1;
	double cons2;
	double cons3;
	double cons4;

	//(mom)->operator()(i,i1) = c5*c2*((mom)->operator()(i,i1)) + (c5*(c3)+q)*F(i,i1) + (c5*(c4)+r)*R(i,i1);
	cons1 = (*obj).getc5()*(*obj).getc2();
	cons2 = (*obj).getc5()*(*obj).getc3()+(*obj).getq();
	cons3 = (*obj).getc5()*(*obj).getc4()+(*obj).getr();
	cons4 = (*obj).getc1();


	// vector<matrix<double> > savef1forces;
	// vector<matrix<double> > savef2forces;
	// vector<matrix<double> > savef3forces;
	// vector<matrix<double> > savef4forces;
	// vector<matrix<double> > savef5forces;
	// vector<matrix<double> > savef6forces;
	// vector<matrix<double> > savef7forces;

	// vector<matrix<double> > saveftemp1forces;
	// vector<matrix<double> > saveftemp2forces;
	// vector<matrix<double> > saveftemp3forces;

	// vector<matrix<double> > savepositions;
	// vector<vector1<int> > savebound;
	// vector<vector1<double> > saveboundalongs;

	for(i = 0 ; i < runtime ; i++) {
		//cout << i << endl;

		cout << i << endl;
		//pausel();
	
		//cout << (*obj).avmom() << endl;
	if(i%25==0) {
		//delete froyo1,froyo2,froyo3,froyo4,froyo5,froyo6;
		// cout << "updated after: " << i << endl;
		// state = obj->getdat();

		delete froyo1,froyo2,froyo3,froyo4,froyo5,froyo6;
		// cout << "updated after: " << i << endl;
		// state = obj->getdat();
		froyo1 = obj->calculatepairs(boxes,pai,3.5);
		froyo2 = obj->calculatepairs(boxes,pbi,3.5);
		froyo3 = obj->calculatepairs(boxes,pci,3.5);
		froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
		froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
		froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);


		cudaFree(d1_list1);cudaFree(d1_list2);cudaFree(d1_list3);cudaFree(d1_list4);
		

		cudaFree(d2_list1);cudaFree(d2_list2);cudaFree(d2_list3);cudaFree(d2_list4);
		

		cudaFree(d3_list1);cudaFree(d3_list2);cudaFree(d3_list3);cudaFree(d3_list4);
		

		cudaFree(d4_list1);cudaFree(d4_list2);cudaFree(d4_list3);cudaFree(d4_list4);
		

		cudaFree(d5_list1);cudaFree(d5_list2);cudaFree(d5_list3);cudaFree(d5_list4);
	
		cudaFree(d6_list1);cudaFree(d6_list2);cudaFree(d6_list3);cudaFree(d6_list4);

		cudaFree(d_indices1);

		cudaFree(d_indices2);

		cudaFree(d_close);



		this->resetindices(d_p_indices,totalN);



		construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp,false);


		cout << "pair" << endl;

		pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa
		pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2); //fbb
		pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3); //fcc
		pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4); //fab
		pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5); //fac
		pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc


		
		// froyo1 = obj->calculatepairs(boxes,pai,3.5);
		// froyo2 = obj->calculatepairs(boxes,pbi,3.5);
		// froyo3 = obj->calculatepairs(boxes,pci,3.5);
		// froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
		// froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
		// froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);

	}





	//cout << "pairs" << endl;
	// cout << "pairings" << endl;


	// matrix<double> ftemp2(totalN,dimension),ftemp3(totalN,dimension);
	// //matrix<double> angforces1(nc,dimension-1),angforces2(nc,dimension-1),angforces3(nc,dimension-1);
	// // cout << "matrices initialized" << endl;

	// matrix<double> F1((*obj).calculateforces(*froyo1,*faa)); //calculate the forces using the pairs as an input

	// matrix<double> F2((*obj).calculateforces(*froyo2,*fbb)); //calculate the forces using the pairs as an input

	// matrix<double> ftemp1((*obj).calculateforces(*froyo3,*fcc)); //calculate the forces using the pairs as an input
	
	// matrix<double> F3((*obj).calculateforces(*froyo4,*fab)); //calculate the forces using the pairs as an input

	// this->ForcesDueToPositionPL(*froyo5,ftemp2); //calculate the forces using the pairs as an input

	// this->ForcesDueToPositionPL(*froyo6,ftemp3); //calculate the forces using the pairs as an input

	// this->CalculateBindings(*froyo5,*froyo6);

	// matrix<double> F4 = this->BindingForces();

	// matrix<double> F5 = this->PositionForcesDueToAngles();

	//print_device_float2(d_particles,totalN);

	double *d1_forces1x;
	double *d1_forces2x;
	double *d1_forces1y;
	double *d1_forces2y;
	cudaMalloc((void**)&d1_forces1x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces1y,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2y,th1*sizeof(double));
	calculateforces2D(d1_list1,d1_list2,d_particles, d1_forces1x,d1_forces1y,d1_forces2x,d1_forces2y, faa_gpu ,th1, l,true);

	cout << "force1" << endl;
	arracychck(d1_forces1x,th1);
	arracychck(d1_forces1y,th1); 

	double *d2_forces1x;
	double *d2_forces2x;
	double *d2_forces1y;
	double *d2_forces2y;
	cudaMalloc((void**)&d2_forces1x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces1y,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2y,th2*sizeof(double));
	calculateforces2D(d2_list1,d2_list2,d_particles, d2_forces1x,d2_forces1y,d2_forces2x,d2_forces2y, fbb_gpu ,th2, l,true);	

	cout << "force2" << endl;
	arracychck(d1_forces2x,th2);
	arracychck(d1_forces2y,th2); 

	double *d3_forces1x;
	double *d3_forces2x;
	double *d3_forces1y;
	double *d3_forces2y;
	cudaMalloc((void**)&d3_forces1x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces1y,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2y,th3*sizeof(double));
	calculateforces2D(d3_list1,d3_list2,d_particles, d3_forces1x,d3_forces1y,d3_forces2x,d3_forces2y, fcc_gpu ,th3, l,true);	


	cout << "force3" << endl;
	arracychck(d3_forces1x,th3);
	arracychck(d3_forces1y,th3); 

	double *d4_forces1x;
	double *d4_forces2x;
	double *d4_forces1y;
	double *d4_forces2y;
	cudaMalloc((void**)&d4_forces1x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces1y,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2y,th4*sizeof(double));
	calculateforces2D(d4_list1,d4_list2,d_particles, d4_forces1x,d4_forces1y,d4_forces2x,d4_forces2y, fab_gpu ,th4, l,true);

	cout << "force4" << endl;
	arracychck(d4_forces1x,th4);
	arracychck(d4_forces1y,th4); 

	double *d5_forces1x;
	double *d5_forces2x;
	double *d5_forces1y;
	double *d5_forces2y;
	cudaMalloc((void**)&d5_forces1x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces1y,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2y,th5*sizeof(double));
	calculateforces2D(d5_list1,d5_list2,d_particles, d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y, fbc_gpu ,th5, l,true);

	cout << "force5" << endl;
	arracychck(d5_forces1x,th5);
	arracychck(d5_forces1y,th5); 	
	// cout << "d5" << endl;
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y,d_particles,th5,totalN);
	// pausel();	

	double *d6_forces1x;
	double *d6_forces2x;
	double *d6_forces1y;
	double *d6_forces2y;
	cudaMalloc((void**)&d6_forces1x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces1y,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2y,th6*sizeof(double));
	calculateforces2D(d6_list1,d6_list2,d_particles, d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y, fac_gpu ,th6, l,true);	

	cout << "force6" << endl;
	arracychck(d6_forces1x,th6);
	arracychck(d6_forces1y,th6); 	
	// cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();

	// matrix<double> F6((*obj).calculateforces(*bondpairs,*bindm));
	callCalculateUnbindingsGPU(d_particles,d_bound,d_boundalong,d_changestate);



	cout << "unbindings calculated" << endl;




	callCalculateBindingsGPU(d5_list1,d5_list2,d6_list1,d6_list2,d_particles, d_bound, d_boundalong,d_changestate,th5 ,th6 );


	cout << "bindings calculated" << endl;



	int *d7_list1,*d7_list2,*d7_list3;
	double *d7_forces1x;
	double *d7_forces1y;
	double *d7_forces2x;
	double *d7_forces2y;
	double *d7_forces3x;
	double *d7_forces3y;
	int th7;
	BindingForcesGPU(d_particles, d_bound, d_boundalong, d7_list1,d7_list2,d7_list3, d7_forces1x, d7_forces1y, d7_forces2x,d7_forces2y,d7_forces3x, d7_forces3y, bindp_gpu, th7);


	cout << "binding forces calculated" << endl;
	arracychck(d7_forces1x,th7);
	arracychck(d7_forces1y,th7); 

	double *d8_forces1x;
	double *d8_forces2x;
	double *d8_forces1y;
	double *d8_forces2y;
	cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));
	calculateforces2D(d8_list1,d8_list2,d_particles, d8_forces1x,d8_forces1y,d8_forces2x,d8_forces2y, bindm_gpu ,th8, l,true);

	cout << "forces8" << endl;
	arracychck(d8_forces1x,th8);
	arracychck(d8_forces1y,th8); 

	resetchangestate(d_changestate);
	// matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets,*bendp));

//	int th9;


	double *d9_forces1x;
	double *d9_forces2x;
	double *d9_forces1y;
	double *d9_forces2y;
	double *d9_forces3x;
	double *d9_forces3y;
	cudaMalloc((void**)&d9_forces1x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces1y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3y,th9*sizeof(double));
	BendingForcesGPU(d_particles, d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces1y,d9_forces2x,d9_forces2y, d9_forces3x, d9_forces3y,bendp_gpu,th9);

	cout << "bending forces" << endl;
	arracychck(d9_forces1x,th9);
	arracychck(d9_forces1y,th9); 
	//	matrix<double> F = ftemp1+ftemp2+ftemp3+F1+F2+F3+F4+F5+F6+F7;//+F4+F5;

		// matrix<double> R(totalN,dimension);
		// for(int i1 = 0 ; i1 < totalN ; i1++) {
		// 	for(int j = 0 ; j < dimension ; j++) {
		// 		R(i1,j) = (3.464101615 * ((double) rand() / (RAND_MAX)) - 1.732050808);
		// 	}
		// }

	int *d10_list1;
	double *d10_forces1x;
	double *d10_forces1y;	
	int th10;

	PositionForcesDueToAnglesGPU(d_particles, d_bound, d_boundalong, d10_list1, d10_forces1x, d10_forces1y,th10);
	

	cout << "pos forces" << endl;
	arracychck(d10_forces1x,th9);
	arracychck(d10_forces1y,th9); 

	cout << "all forces calculated" << endl;




	double *d_totalforcex1;
	double *d_totalforcey1;

	cudaMalloc((void**)&d_totalforcex1,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey1,totalN*sizeof(double));

	resetforce(d_totalforcex1);
	resetforce(d_totalforcey1);

	cout << "reset" << endl;

	// print_device_array(d_totalforcex,totalN);


	//ReduceForces(d1_list1,d1_list2,d1_list3,d1_list4,d1_forces1x,d1_forces2x,d1_forces1y,d1_forces2y,d_totalforcex,d_totalforcey,th1);
	ReduceForces(d1_list1,d1_list2,d1_list3,d1_list4,d1_forces1x,d1_forces2x,d1_forces1y,d1_forces2y,d_totalforcex1,d_totalforcey1,th1);	
	 cout << "d1" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	

	double *d_totalforcex2;
	double *d_totalforcey2;

	cudaMalloc((void**)&d_totalforcex2,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey2,totalN*sizeof(double));

	resetforce(d_totalforcex2);
	resetforce(d_totalforcey2);	 
	//ReduceForces(d2_list1,d2_list2,d2_list3,d2_list4,d2_forces1x,d2_forces2x,d2_forces1y,d2_forces2y,d_totalforcex,d_totalforcey,th2);
	ReduceForces(d2_list1,d2_list2,d2_list3,d2_list4,d2_forces1x,d2_forces2x,d2_forces1y,d2_forces2y,d_totalforcex2,d_totalforcey2,th2);
	cout << "d2" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex3;
	double *d_totalforcey3;

	cudaMalloc((void**)&d_totalforcex3,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey3,totalN*sizeof(double));

	resetforce(d_totalforcex3);
	resetforce(d_totalforcey3);	 	
	//ReduceForces(d3_list1,d3_list2,d3_list3,d3_list4,d3_forces1x,d3_forces2x,d3_forces1y,d3_forces2y,d_totalforcex,d_totalforcey,th3);	 	
	ReduceForces(d3_list1,d3_list2,d3_list3,d3_list4,d3_forces1x,d3_forces2x,d3_forces1y,d3_forces2y,d_totalforcex3,d_totalforcey3,th3);
	 cout << "d3" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex4;
	double *d_totalforcey4;

	cudaMalloc((void**)&d_totalforcex4,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey4,totalN*sizeof(double));

	resetforce(d_totalforcex4);
	resetforce(d_totalforcey4);	 		
	//ReduceForces(d4_list1,d4_list2,d4_list3,d4_list4,d4_forces1x,d4_forces2x,d4_forces1y,d4_forces2y,d_totalforcex,d_totalforcey,th4);
	ReduceForces(d4_list1,d4_list2,d4_list3,d4_list4,d4_forces1x,d4_forces2x,d4_forces1y,d4_forces2y,d_totalforcex4,d_totalforcey4,th4);
	 cout << "d4" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex5;
	double *d_totalforcey5;

	cudaMalloc((void**)&d_totalforcex5,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey5,totalN*sizeof(double));

	resetforce(d_totalforcex5);
	resetforce(d_totalforcey5);	 	 
	//ReduceForces(d5_list1,d5_list2,d5_list3,d5_list4,d5_forces1x,d5_forces2x,d5_forces1y,d5_forces2y,d_totalforcex,d_totalforcey,th5);
	ReduceForces(d5_list1,d5_list2,d5_list3,d5_list4,d5_forces1x,d5_forces2x,d5_forces1y,d5_forces2y,d_totalforcex5,d_totalforcey5,th5);	

	 cout << "d5" << endl;	
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d_particles,th5,totalN);
	// pausel();
	double *d_totalforcex6;
	double *d_totalforcey6;

	cudaMalloc((void**)&d_totalforcex6,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey6,totalN*sizeof(double));

	resetforce(d_totalforcex6);
	resetforce(d_totalforcey6);	 	
	//ReduceForces(d6_list1,d6_list2,d6_list3,d6_list4,d6_forces1x,d6_forces2x,d6_forces1y,d6_forces2y,d_totalforcex,d_totalforcey,th6);	 	
	ReduceForces(d6_list1,d6_list2,d6_list3,d6_list4,d6_forces1x,d6_forces2x,d6_forces1y,d6_forces2y,d_totalforcex6,d_totalforcey6,th6);		
	 cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();
	double *d_totalforcex7;
	double *d_totalforcey7;

	cudaMalloc((void**)&d_totalforcex7,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey7,totalN*sizeof(double));

	resetforce(d_totalforcex7);
	resetforce(d_totalforcey7);	 
	//ReduceForces3(d7_list1,d7_list2,d7_list3,d7_forces1x,d7_forces2x,d7_forces3x,d7_forces1y,d7_forces2y,d7_forces3y,d_totalforcex,d_totalforcey,th7);	 
	ReduceForces3(d7_list1,d7_list2,d7_list3,d7_forces1x,d7_forces2x,d7_forces3x,d7_forces1y,d7_forces2y,d7_forces3y,d_totalforcex7,d_totalforcey7,th7);	
	 cout << "d7" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex8;
	double *d_totalforcey8;

	cudaMalloc((void**)&d_totalforcex8,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey8,totalN*sizeof(double));
	

	
	resetforce(d_totalforcex8);
	resetforce(d_totalforcey8);	 
	//ReduceForces(d8_list1,d8_list2,d8_forces1x,d8_forces2x,d8_forces1y,d8_forces2y,d_totalforcex,d_totalforcey,th8);	 
	ReduceForces(d8_list1,d8_list2,d8_forces1x,d8_forces2x,d8_forces1y,d8_forces2y,d_totalforcex8,d_totalforcey8,th8);
	 cout << "d8" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex9;
	double *d_totalforcey9;

	cudaMalloc((void**)&d_totalforcex9,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey9,totalN*sizeof(double));

	resetforce(d_totalforcex9);
	resetforce(d_totalforcey9);	 
	//ReduceForces3(d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces2x,d9_forces3x,d9_forces1y,d9_forces2y,d9_forces3y,d_totalforcex,d_totalforcey,th9);		 
	ReduceForces3(d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces2x,d9_forces3x,d9_forces1y,d9_forces2y,d9_forces3y,d_totalforcex9,d_totalforcey9,th9);		
	
	double *d_totalforcex10;
	double *d_totalforcey10;

	cudaMalloc((void**)&d_totalforcex10,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey10,totalN*sizeof(double));

	resetforce(d_totalforcex10);
	resetforce(d_totalforcey10);	 
	 cout << "d9" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double ff = (v0_a+v0_b)/2.;
	//ReduceForcesAndNormalize(d10_list1,d10_forces1x,d10_forces1y,d_totalforcex,d_totalforcey, max_s, ff, th10);
	ReduceForcesAndNormalize(d10_list1,d10_forces1x,d10_forces1y,d_totalforcex10,d_totalforcey10, max_s, ff, th10);

	cout << "reduction" << endl;
	// cout << "d10" << endl;
	// print_device_array(d_totalforcex,totalN);
	// pausel();

	matrix<double> ftemp2(totalN,dimension),ftemp3(totalN,dimension);
	//matrix<double> angforces1(nc,dimension-1),angforces2(nc,dimension-1),angforces3(nc,dimension-1);
	// cout << "matrices initialized" << endl;

	matrix<double> F1((*obj).calculateforces(*froyo1,*faa)); //calculate the forces using the pairs as an input

	matrix<double> F2((*obj).calculateforces(*froyo2,*fbb)); //calculate the forces using the pairs as an input

	matrix<double> ftemp1((*obj).calculateforces(*froyo3,*fcc)); //calculate the forces using the pairs as an input
	
	matrix<double> F3((*obj).calculateforces(*froyo4,*fab)); //calculate the forces using the pairs as an input

	this->ForcesDueToPositionPL(*froyo5,ftemp2); //calculate the forces using the pairs as an input

	this->ForcesDueToPositionPL(*froyo6,ftemp3); //calculate the forces using the pairs as an input

	this->CalculateBindings(*froyo5,*froyo6);

	matrix<double> F4 = this->BindingForces();

	matrix<double> F5 = this->PositionForcesDueToAngles();


//	cout << "active forces" << endl;
	//cout << "pos forces" << endl;
	 matrix<double> F6((*obj).calculateforces(*bondpairs,*bindm));

	//cout << "bond forces" << endl;
	matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets,*bendp));

	//cout << "after check matrix" << endl;

	matrix<double> F = ftemp1+ftemp2+ftemp3+F1+F2+F3+F4+F5+F6+F7;

	cout << l << endl;
	for(int j1 = 0 ; j1 < na+nb ; j1++ ) {
		if(bound[j1]>0) cout << j1 << ",";
	}
	cout << endl;
	print_device_array_indices(d_bound,na+nb);
	pausel();




	cout << F1 << endl;
	print_device_array_weave(d_totalforcex1,d_totalforcey1,totalN);
	cout << th1 << endl;
	cout << "faa" << endl;
	pausel();



	cout << F2 << endl;
	print_device_array_weave(d_totalforcex2,d_totalforcey2,totalN);
	cout << "fbb" << endl;
	pausel();


	cout << ftemp1 << endl;
	print_device_array_weave(d_totalforcex3,d_totalforcey3,totalN);
	cout << "fcc" << endl;
	pausel();		

	cout << F3 << endl;
	print_device_array_weave(d_totalforcex4,d_totalforcey4,totalN);
	cout << "fab" << endl;
	pausel();

	cout << ftemp2 << endl;
	print_device_array_weave(d_totalforcex5,d_totalforcey5,totalN);
	cout << "fac" << endl;
	pausel();		

	cout << ftemp3 << endl;
	print_device_array_weave(d_totalforcex6,d_totalforcey6,totalN);
	cout << "fbc" << endl;
	pausel();

	cout << F4 << endl;
	print_device_array_weave(d_totalforcex7,d_totalforcey7,totalN);
	cout << "bound to mt" << endl;
	pausel();	

	cout << F5 << endl;
	print_device_array_weave(d_totalforcex10,d_totalforcey10,totalN);
	cout << "position force" << endl;
	pausel();

	cout << F6 << endl;
	print_device_array_weave(d_totalforcex8,d_totalforcey8,totalN);
	cout << "bound within mt" << endl;
	pausel();		


	cout << F7 << endl;
	print_device_array_weave(d_totalforcex9,d_totalforcey9,totalN);
	cout << "bending force" << endl;
	pausel();		







	double *d_R1;
	double *d_R2;

	// cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));	
	cudaMalloc((void**)&d_R1,totalN*sizeof(double));
	cudaMalloc((void**)&d_R2,totalN*sizeof(double));

	setstaterandom(d_R1,1.732050808,totalN);
 	setstaterandom(d_R2,1.732050808,totalN);

 	print_device_array(d_R1,totalN);
 	print_device_array(d_R2,totalN);

 	pausel();

 	advmom2D(d_momenta, d_totalforcex, d_totalforcey, d_R1, d_R2, cons1,cons2,cons3,totalN);
 	advpos2D(d_particles, d_momenta, cons4, totalN);

 	applypbc2D(d_particles,d_momenta,l,is_periodic,totalN);

 	cout << "updated" << endl;

	cudaFree(d1_forces1x);
	cudaFree(d1_forces2x);
	cudaFree(d1_forces1y);
	cudaFree(d1_forces2y);
	cudaFree(d2_forces1x);
	cudaFree(d2_forces2x);
	cudaFree(d2_forces1y);
	cudaFree(d2_forces2y);
	cudaFree(d3_forces1x);
	cudaFree(d3_forces2x);
	cudaFree(d3_forces1y);
	cudaFree(d3_forces2y);
	cudaFree(d4_forces1x);
	cudaFree(d4_forces2x);
	cudaFree(d4_forces1y);
	cudaFree(d4_forces2y);
	cudaFree(d5_forces1x);
	cudaFree(d5_forces2x);
	cudaFree(d5_forces1y);
	cudaFree(d5_forces2y);
	cudaFree(d6_forces1x);
	cudaFree(d6_forces2x);
	cudaFree(d6_forces1y);
	cudaFree(d6_forces2y);
 	cudaFree(d7_list1);
 	cudaFree(d7_list2);
 	cudaFree(d7_list3);
	cudaFree(d7_forces1x);
	cudaFree(d7_forces2x);
	cudaFree(d7_forces3x);
	cudaFree(d7_forces1y);
	cudaFree(d7_forces2y);
	cudaFree(d7_forces3y);
	cudaFree(d8_forces1x);
	cudaFree(d8_forces2x);
	cudaFree(d8_forces1y);
	cudaFree(d8_forces2y);
	cudaFree(d9_forces1x);
	cudaFree(d9_forces2x);
	cudaFree(d9_forces3x);
	cudaFree(d9_forces1y);
	cudaFree(d9_forces2y);
	cudaFree(d9_forces3y);
	cudaFree(d10_list1);
	cudaFree(d10_forces1x);
	cudaFree(d10_forces1y);

	cout << "freed" << endl;


	}
}
*/

template <typename Fun>
void Microtubule::runGPUSV(int runtime, int every, Fun func)
{
//	pausel();
	//int ccc;
	int totalN = obj->getN();

	//num is the number of boxes per length

	int ncells = num*num;

	WCApotentialGPU faa_gpu(2.,1.,2.);
	WCApotentialGPU fab_gpu(1.,1.,0.);
	WCApotentialGPU fac_gpu(1.,1.,0.);
	WCApotentialGPU fbb_gpu(2.,1.,2.);
	WCApotentialGPU fbc_gpu(1.,1.,0.);
	WCApotentialGPU fcc_gpu(1.,1.,0.);
	HarmonicPotentialGPU bindp_gpu(100.,1.); 
	FENEPotentialGPU bindm_gpu(50.,1.5); 
	BendingPotentialGPU bendp_gpu(100.,0.);

//we now have the count of each cell list
	int nbpairs = 5*ncells;
	int nperl = num;

	int *cells1 = new int [nbpairs];
	int *cells2 = new int [nbpairs];


	int itery = 0;
	for(int i1 = 0 ; i1 < num ; i1++) {
		for(int i2 = 0 ; i2 < num ; i2++ ) {


			int b1 =  i1*nperl+i2;

			int i3 = i1+0;
			int j3 = i2+0;

			int i4 = i1+1;
			int j4 = i2+0;

			int i5 = i1-1;
			int j5 = i2+1;

			int i6 = i1+0;
			int j6 = i2+1;

			int i7 = i1+1;
			int j7 = i2+1;

			prdshft(i3,nperl);
			prdshft(j3,nperl);

			prdshft(i4,nperl);
			prdshft(j4,nperl);

			prdshft(i5,nperl);
			prdshft(j5,nperl);
			
			prdshft(i6,nperl);
			prdshft(j6,nperl);
			
			prdshft(i7,nperl);
			prdshft(j7,nperl);		

			cells1[itery] =  b1;
			cells2[itery] =  i3*nperl+j3;

			itery++;

			cells1[itery] =  b1;
			cells2[itery] =  i4*nperl+j4;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i5*nperl+j5;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i6*nperl+j6;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i7*nperl+j7;

			itery++;


		}
	}
	int size4 = nbpairs*sizeof(int);

	int *d_cells1;
	int *d_cells2;

	cudaMalloc((void**)&d_cells1,size4);

	cudaMalloc((void**)&d_cells2,size4);

	cudaMemcpy(d_cells1,cells1,size4,cudaMemcpyHostToDevice);
	cudaMemcpy(d_cells2,cells2,size4,cudaMemcpyHostToDevice);
	//int sibdiv = floor(ll/4.0);
	// print_device_array(d_cells1,nbpairs);
	// print_device_array(d_cells2,nbpairs);



	// matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);
	
	float2 *particles = new float2 [totalN];
	float2 *momenta = new float2 [totalN];
	int *p_indices = new int [totalN];

	for(int i = 0 ; i < totalN ; i++)
	p_indices[i]=i;

	float2 *d_particles;
	float2 *d_momenta;
	int *d_p_indices;

	int *d_bound;
	double *d_boundalong;
	int *d_changestate;

	cudaMalloc((void**)&d_bound,(na+nb)*sizeof(int));
	cudaMalloc((void**)&d_boundalong,(na+nb)*sizeof(double));
	cudaMalloc((void**)&d_changestate,(na+nb)*sizeof(int));

	cudaMemset(d_bound,0,(na+nb)*sizeof(int));
	cudaMemset(d_boundalong,0.,(na+nb)*sizeof(double));
	cudaMemset(d_changestate,0,(na+nb)*sizeof(int));

	double *d_totalforcex;
	double *d_totalforcey;

	cudaMalloc((void**)&d_totalforcex,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey,totalN*sizeof(double));


	matrix<double> state(obj->getdat());


	for(int i = 0  ; i < totalN ; i++) {

	float2 c;
	c.x=state(i,0);
	c.y=state(i,1);

	(particles)[i]=c;

	float2 d;

	d.x = 0.;
	d.y = 0.;

	(momenta)[i]=d;
	}


	int size =  totalN*sizeof(float2);
	int size2 = totalN*sizeof(int);


	cudaMalloc((void**)&d_particles,size);
	cudaMalloc((void**)&d_momenta,size);
	cudaMalloc((void**)&d_p_indices,size2);

	cudaMemcpy(d_particles,particles,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_momenta,momenta,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_p_indices,p_indices,size2,cudaMemcpyHostToDevice);


	// matrix<int> *froyo1 = obj->calculatepairs(boxes,pai,3.5);
	// matrix<int> *froyo2 = obj->calculatepairs(boxes,pbi,3.5);
	// matrix<int> *froyo3 = obj->calculatepairs(boxes,pci,3.5);
	// matrix<int> *froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
	// matrix<int> *froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
	// matrix<int> *froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);
	
	int *d_indices1;
	int *d_indices2;
	double *d_close;


	int tpp;
	


	construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp);



	less_than_condition_AND cond1(SQR(3.5),0,na);
	less_than_condition_AND cond2(SQR(3.5),na,na+nb);
	less_than_condition_AND cond3(SQR(3.5),na+nb,na+nb+nc);
	less_than_condition_NAND cond4(SQR(3.5),0,na,na,na+nb);
	less_than_condition_NAND cond5(SQR(3.5),0,na,na+nb,na+nb+nc);
	less_than_condition_NAND cond6(SQR(3.5),na,na+nb,na+nb,na+nb+nc);


	int th1;
	int *d1_list1,*d1_list2,*d1_list3,*d1_list4;
	pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa

	int th2;
	int *d2_list1,*d2_list2,*d2_list3,*d2_list4;
	pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2);	//fbb

	int th3;
	int *d3_list1,*d3_list2,*d3_list3,*d3_list4;
	pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3);	//fcc

	int th4;
	int *d4_list1,*d4_list2,*d4_list3,*d4_list4;
	pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4);	//fab

	int th5;
	int *d5_list1,*d5_list2,*d5_list3,*d5_list4;
	pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5);	//fac

	int th6;
	int *d6_list1,*d6_list2,*d6_list3,*d6_list4;
	pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc			
	//matrix<double> state(obj->getdat()); //the state of the system

	int th8 = (*bondpairs).getNsafe();
	int th9 = (*bendtriplets).getNsafe();


	int *d8_list1,*d8_list2,*d8_list3,*d8_list4;
	cudaMalloc((void**)&d8_list1,th8*sizeof(int));
	cudaMalloc((void**)&d8_list2,th8*sizeof(int));
	cudaMalloc((void**)&d8_list3,th8*sizeof(int));
	cudaMalloc((void**)&d8_list4,th8*sizeof(int));



	int *h8_list1 = new int [th8];
	int *h8_list2 = new int [th8];

	for(int i = 0 ; i < th8 ; i++ ) {
		h8_list1[i] = (*bondpairs)(i,0);
		h8_list2[i] = (*bondpairs)(i,1);
	}

	cudaMemcpy(d8_list1,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list2,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list3,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list4,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);





	int *d9_list1,*d9_list2,*d9_list3,*d9_list4,*d9_list5,*d9_list6;
	cudaMalloc((void**)&d9_list1,th9*sizeof(int));
	cudaMalloc((void**)&d9_list2,th9*sizeof(int));
	cudaMalloc((void**)&d9_list3,th9*sizeof(int));
	cudaMalloc((void**)&d9_list4,th9*sizeof(int));
	cudaMalloc((void**)&d9_list5,th9*sizeof(int));
	cudaMalloc((void**)&d9_list6,th9*sizeof(int));

	int *h9_list1 = new int [th9];
	int *h9_list2 = new int [th9];
	int *h9_list3 = new int [th9];

	for(int i = 0 ; i < th9 ; i++ ) {
		h9_list1[i] = (*bendtriplets)(i,0);
		h9_list2[i] = (*bendtriplets)(i,1);
		h9_list3[i] = (*bendtriplets)(i,2);
	}	

	cudaMemcpy(d9_list1,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list2,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list3,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list4,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list5,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list6,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);


	int i;

	double cons1;
	double cons2;
	double cons3;
	double cons4;

	//(mom)->operator()(i,i1) = c5*c2*((mom)->operator()(i,i1)) + (c5*(c3)+q)*F(i,i1) + (c5*(c4)+r)*R(i,i1);
	cons1 = (*obj).getc5()*(*obj).getc2();
	cons2 = (*obj).getc5()*(*obj).getc3()+(*obj).getq();
	cons3 = (*obj).getc5()*(*obj).getc4()+(*obj).getr();
	cons4 = (*obj).getc1();

	double d_dt = (*obj).getdt();
	double d_m = (*obj).getm();
	double d_kT = (*obj).getkT();


	// vector<matrix<double> > savef1forces;
	// vector<matrix<double> > savef2forces;
	// vector<matrix<double> > savef3forces;
	// vector<matrix<double> > savef4forces;
	// vector<matrix<double> > savef5forces;
	// vector<matrix<double> > savef6forces;
	// vector<matrix<double> > savef7forces;

	// vector<matrix<double> > saveftemp1forces;
	// vector<matrix<double> > saveftemp2forces;
	// vector<matrix<double> > saveftemp3forces;

	// vector<matrix<double> > savepositions;
	// vector<vector1<int> > savebound;
	// vector<vector1<double> > saveboundalongs;

	for(i = 0 ; i < runtime ; i++) {
		//cout << i << endl;

		cout << i << endl;
	
		//cout << (*obj).avmom() << endl;
	if(i%25==0) {
		//delete froyo1,froyo2,froyo3,froyo4,froyo5,froyo6;
		// cout << "updated after: " << i << endl;
		// state = obj->getdat();


		cudaFree(d1_list1);cudaFree(d1_list2);cudaFree(d1_list3);cudaFree(d1_list4);
		

		cudaFree(d2_list1);cudaFree(d2_list2);cudaFree(d2_list3);cudaFree(d2_list4);
		

		cudaFree(d3_list1);cudaFree(d3_list2);cudaFree(d3_list3);cudaFree(d3_list4);
		

		cudaFree(d4_list1);cudaFree(d4_list2);cudaFree(d4_list3);cudaFree(d4_list4);
		

		cudaFree(d5_list1);cudaFree(d5_list2);cudaFree(d5_list3);cudaFree(d5_list4);
	
		cudaFree(d6_list1);cudaFree(d6_list2);cudaFree(d6_list3);cudaFree(d6_list4);

		cudaFree(d_indices1);

		cudaFree(d_indices2);

		cudaFree(d_close);

		this->resetindices(d_p_indices,totalN);


		construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp,false);

		// cout << tpp << endl;
		// cout << "pair" << endl;

		pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa
		pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2); //fbb
		pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3); //fcc
		pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4); //fab
		pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5); //fac
		pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc

	}




	double *d1_forces1x;
	double *d1_forces2x;
	double *d1_forces1y;
	double *d1_forces2y;
	cudaMalloc((void**)&d1_forces1x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces1y,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2y,th1*sizeof(double));
	calculateforces2D(d1_list1,d1_list2,d_particles, d1_forces1x,d1_forces1y,d1_forces2x,d1_forces2y, faa_gpu ,th1, l,true);

	// arracychck(d1_forces1x,th1);
	// arracychck(d1_forces1y,th1);
	// cout << "force1" << endl;

	double *d2_forces1x;
	double *d2_forces2x;
	double *d2_forces1y;
	double *d2_forces2y;
	cudaMalloc((void**)&d2_forces1x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces1y,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2y,th2*sizeof(double));
	calculateforces2D(d2_list1,d2_list2,d_particles, d2_forces1x,d2_forces1y,d2_forces2x,d2_forces2y, fbb_gpu ,th2, l,true);	



	// arracychck(d2_forces1x,th2);
	// arracychck(d2_forces1y,th2);

	// cout << "force2" << endl;
	double *d3_forces1x;
	double *d3_forces2x;
	double *d3_forces1y;
	double *d3_forces2y;
	cudaMalloc((void**)&d3_forces1x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces1y,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2y,th3*sizeof(double));
	calculateforces2D(d3_list1,d3_list2,d_particles, d3_forces1x,d3_forces1y,d3_forces2x,d3_forces2y, fcc_gpu ,th3, l,true);	

	// arracychck(d3_forces1x,th3);
	// arracychck(d3_forces1y,th3);

	// cout << "force3" << endl;
	double *d4_forces1x;
	double *d4_forces2x;
	double *d4_forces1y;
	double *d4_forces2y;
	cudaMalloc((void**)&d4_forces1x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces1y,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2y,th4*sizeof(double));
	calculateforces2D(d4_list1,d4_list2,d_particles, d4_forces1x,d4_forces1y,d4_forces2x,d4_forces2y, fab_gpu ,th4, l,true);
	
	// arracychck(d4_forces1x,th4);
	// arracychck(d4_forces1y,th4);
	// cout << "force4" << endl;

	double *d5_forces1x;
	double *d5_forces2x;
	double *d5_forces1y;
	double *d5_forces2y;
	cudaMalloc((void**)&d5_forces1x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces1y,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2y,th5*sizeof(double));
	calculateforces2D(d5_list1,d5_list2,d_particles, d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y, fbc_gpu ,th5, l,true);
	// arracychck(d5_forces1x,th5);
	// arracychck(d5_forces1y,th5);
	// cout << "force5" << endl;
	// cout << "d5" << endl;
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y,d_particles,th5,totalN);
	// pausel();	

	double *d6_forces1x;
	double *d6_forces2x;
	double *d6_forces1y;
	double *d6_forces2y;
	cudaMalloc((void**)&d6_forces1x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces1y,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2y,th6*sizeof(double));
	calculateforces2D(d6_list1,d6_list2,d_particles, d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y, fac_gpu ,th6, l,true);	
	// arracychck(d6_forces1x,th6);
	// arracychck(d6_forces1y,th6);
	// cout << "force6" << endl;
	// cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();

	// matrix<double> F6((*obj).calculateforces(*bondpairs,*bindm));
	callCalculateUnbindingsGPU(d_particles,d_bound,d_boundalong,d_changestate);



	// cout << "unbindings calculated" << endl;


	callCalculateBindingsGPU(d5_list1,d5_list2,d6_list1,d6_list2,d_particles, d_bound, d_boundalong,d_changestate,th5 ,th6 );


	// cout << "bindings calculated" << endl;


	int *d7_list1,*d7_list2,*d7_list3;
	double *d7_forces1x;
	double *d7_forces1y;
	double *d7_forces2x;
	double *d7_forces2y;
	double *d7_forces3x;
	double *d7_forces3y;
	int th7;
	BindingForcesGPU(d_particles, d_bound, d_boundalong, d7_list1,d7_list2,d7_list3, d7_forces1x, d7_forces1y, d7_forces2x,d7_forces2y,d7_forces3x, d7_forces3y, bindp_gpu, th7);

	// arracychck(d7_forces1x,th7);
	// arracychck(d7_forces1y,th7);
	// cout << "binding forces calculated" << endl;


	double *d8_forces1x;
	double *d8_forces2x;
	double *d8_forces1y;
	double *d8_forces2y;
	cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));
	calculateforces2D(d8_list1,d8_list2,d_particles, d8_forces1x,d8_forces1y,d8_forces2x,d8_forces2y, bindm_gpu ,th8, l,true);
	
	// arracychck(d8_forces1x,th8);
	// arracychck(d8_forces1y,th8);

	// cout << "force 8" << endl;
	resetchangestate(d_changestate);
	// matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets,*bendp));

//	int th9;


	double *d9_forces1x;
	double *d9_forces2x;
	double *d9_forces1y;
	double *d9_forces2y;
	double *d9_forces3x;
	double *d9_forces3y;
	// cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));	
	cudaMalloc((void**)&d9_forces1x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces1y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3y,th9*sizeof(double));
	BendingForcesGPU(d_particles, d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces1y,d9_forces2x,d9_forces2y, d9_forces3x, d9_forces3y,bendp_gpu,th9);

	// print_device_array_weave(d9_forces1x,d9_forces1y,th9);
	// print_device_array_weave(d9_forces2x,d9_forces2y,th9);
	// print_device_array_weave(d9_forces3x,d9_forces3y,th9);
	// cout << "force 9" << endl;	
	// arracychck(d9_forces1x,th9);
	// arracychck(d9_forces1y,th9);

	// print_device_float2(d_particles,totalN);
	// print_device_array(d9_list1,th9);
	// print_device_array(d9_list2,th9);
	// print_device_array(d9_list3,th9);

	// cout << "force9" << endl;

	//	matrix<double> F = ftemp1+ftemp2+ftemp3+F1+F2+F3+F4+F5+F6+F7;//+F4+F5;

		// matrix<double> R(totalN,dimension);
		// for(int i1 = 0 ; i1 < totalN ; i1++) {
		// 	for(int j = 0 ; j < dimension ; j++) {
		// 		R(i1,j) = (3.464101615 * ((double) rand() / (RAND_MAX)) - 1.732050808);
		// 	}
		// }

	int *d10_list1;
	double *d10_forces1x;
	double *d10_forces1y;	
	int th10;

	PositionForcesDueToAnglesGPU(d_particles, d_bound, d_boundalong, d10_list1, d10_forces1x, d10_forces1y,th10);
	

	// arracychck(d10_forces1x,th10);
	// arracychck(d10_forces1y,th10);
	// cout << "force10" << endl;
	// cout << "all forces calculated" << endl;





	resetforce(d_totalforcex);

	resetforce(d_totalforcey);

	//cout << "reset" << endl;

	// print_device_array(d_totalforcex,totalN);


	ReduceForces(d1_list1,d1_list2,d1_list3,d1_list4,d1_forces1x,d1_forces2x,d1_forces1y,d1_forces2y,d_totalforcex,d_totalforcey,th1);
	// cout << "d1" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d2_list1,d2_list2,d2_list3,d2_list4,d2_forces1x,d2_forces2x,d2_forces1y,d2_forces2y,d_totalforcex,d_totalforcey,th2);
	/// cout << "d2" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces(d3_list1,d3_list2,d3_list3,d3_list4,d3_forces1x,d3_forces2x,d3_forces1y,d3_forces2y,d_totalforcex,d_totalforcey,th3);
	// cout << "d3" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d4_list1,d4_list2,d4_list3,d4_list4,d4_forces1x,d4_forces2x,d4_forces1y,d4_forces2y,d_totalforcex,d_totalforcey,th4);
	// cout << "d4" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d5_list1,d5_list2,d5_list3,d5_list4,d5_forces1x,d5_forces2x,d5_forces1y,d5_forces2y,d_totalforcex,d_totalforcey,th5);
	// cout << "d5" << endl;	
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d_particles,th5,totalN);
	// pausel();	
	ReduceForces(d6_list1,d6_list2,d6_list3,d6_list4,d6_forces1x,d6_forces2x,d6_forces1y,d6_forces2y,d_totalforcex,d_totalforcey,th6);	
	// cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();
	ReduceForces3(d7_list1,d7_list2,d7_list3,d7_forces1x,d7_forces2x,d7_forces3x,d7_forces1y,d7_forces2y,d7_forces3y,d_totalforcex,d_totalforcey,th7);	
	// cout << "d7" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces(d8_list1,d8_list2,d8_forces1x,d8_forces2x,d8_forces1y,d8_forces2y,d_totalforcex,d_totalforcey,th8);
	// cout << "d8" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces3(d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces2x,d9_forces3x,d9_forces1y,d9_forces2y,d9_forces3y,d_totalforcex,d_totalforcey,th9);		
	
	// cout << "d9" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double ff = (v0_a+v0_b)/2.;
	ReduceForcesAndNormalize(d10_list1,d10_forces1x,d10_forces1y,d_totalforcex,d_totalforcey, max_s, ff, th10);

	cout << "reduction" << endl;
	// cout << "d10" << endl;
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	if(i>0&&i%every==0) { 
		// for(int j = 0 ; j < na+nb ; j++) {
		// if(bound[j]>0){
		// cout << "printed" << endl;
		// cout << j << endl;
		// cout << bound[j] << endl;
		// cout << bound_along[j] << endl;
		// cout << F5(j,'r') << endl;
		// cout << F4(j,'r') << endl;
		// cout << obj->getcoordinate(j,0) << " " << obj->getcoordinate(j,1) << endl;
		// cout << obj->getcoordinate(200,0) <<  " " << obj->getcoordinate(200,1) << endl;
		// cout << obj->getcoordinate(200+L,0) <<  " " << obj->getcoordinate(200+L,1) << endl;
		// }
		// }
		// cout << F6 << endl;
		// cout << F7 << endl;

		// cout << ftemp2 << endl;
		// cout << ftemp3 << endl;

		stringstream ss2;
		// ss2 <<i/every;
		// string pairlist = "list";

		 stringstream kts;
		 kts << (*obj).getkT();

		 // stringstream epi;
		 // epi << eps;

		 // stringstream epieq;
		 // epieq << eqeps;

		 stringstream len;
		 len << l;

		 string extension =  "_kT="+kts.str()+"_l="+len.str()+".csv";

		stringstream ss;
		ss <<(i/every);
		string filename = "x";
		filename += ss.str();
		filename += extension;

		string momname = "bind";
		momname += ss.str();
		momname += extension;

		string baname = "bind_along";
		baname += ss.str();
		baname += extension;



		ofstream myfile;
		myfile.open(filename.c_str());
		//myfile <<= (*obj).getdat();
		file_print_device_float2(d_particles,totalN,myfile);
		myfile.close();



		}	

	double *d_R1;
	double *d_R2;

	// cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));	
	cudaMalloc((void**)&d_R1,totalN*sizeof(double));
	cudaMalloc((void**)&d_R2,totalN*sizeof(double));

	setstaterandom(d_R1,1.732050808,totalN);
 	setstaterandom(d_R2,1.732050808,totalN);	

 	//advmom2D(d_momenta, d_totalforcex, d_totalforcey, d_R1, d_R2, cons1,cons2,cons3,totalN);
 	advmom2D_spatialdependence(d_momenta,d_particles,d_totalforcex,d_totalforcey,d_R1,d_R2,func, d_dt, d_kT, d_m, totalN);
 	advpos2D(d_particles, d_momenta, cons4, totalN);


 	applypbc2D(d_particles,d_momenta,l,is_periodic,totalN);

 	cout << "updated" << endl;

	cudaFree(d1_forces1x);
	cudaFree(d1_forces2x);
	cudaFree(d1_forces1y);
	cudaFree(d1_forces2y);
	cudaFree(d2_forces1x);
	cudaFree(d2_forces2x);
	cudaFree(d2_forces1y);
	cudaFree(d2_forces2y);
	cudaFree(d3_forces1x);
	cudaFree(d3_forces2x);
	cudaFree(d3_forces1y);
	cudaFree(d3_forces2y);
	cudaFree(d4_forces1x);
	cudaFree(d4_forces2x);
	cudaFree(d4_forces1y);
	cudaFree(d4_forces2y);
	cudaFree(d5_forces1x);
	cudaFree(d5_forces2x);
	cudaFree(d5_forces1y);
	cudaFree(d5_forces2y);
	cudaFree(d6_forces1x);
	cudaFree(d6_forces2x);
	cudaFree(d6_forces1y);
	cudaFree(d6_forces2y);
 	cudaFree(d7_list1);
 	cudaFree(d7_list2);
 	cudaFree(d7_list3);
	cudaFree(d7_forces1x);
	cudaFree(d7_forces2x);
	cudaFree(d7_forces3x);
	cudaFree(d7_forces1y);
	cudaFree(d7_forces2y);
	cudaFree(d7_forces3y);
	cudaFree(d8_forces1x);
	cudaFree(d8_forces2x);
	cudaFree(d8_forces1y);
	cudaFree(d8_forces2y);
	cudaFree(d9_forces1x);
	cudaFree(d9_forces2x);
	cudaFree(d9_forces3x);
	cudaFree(d9_forces1y);
	cudaFree(d9_forces2y);
	cudaFree(d9_forces3y);
	cudaFree(d10_list1);
	cudaFree(d10_forces1x);
	cudaFree(d10_forces1y);

	cout << "freed" << endl;


	}
}


template <typename Fun>
void Microtubule::runGPUPV(int runtime, int every, Fun func)
{
//	pausel();
	//int ccc;
	int totalN = obj->getN();

	//num is the number of boxes per length

	int ncells = num*num;

	WCApotentialGPU faa_gpu(2.,1.,2.);
	WCApotentialGPU fab_gpu(1.,1.,0.);
	WCApotentialGPU fac_gpu(1.,1.,0.);
	WCApotentialGPU fbb_gpu(2.,1.,2.);
	WCApotentialGPU fbc_gpu(1.,1.,0.);
	WCApotentialGPU fcc_gpu(1.,1.,0.);
	HarmonicPotentialGPU bindp_gpu(100.,1.); 
	FENEPotentialGPU bindm_gpu(50.,1.5); 
	BendingPotentialGPU bendp_gpu(100.,0.);

//we now have the count of each cell list
	int nbpairs = 5*ncells;
	int nperl = num;

	int *cells1 = new int [nbpairs];
	int *cells2 = new int [nbpairs];


	int itery = 0;
	for(int i1 = 0 ; i1 < num ; i1++) {
		for(int i2 = 0 ; i2 < num ; i2++ ) {


			int b1 =  i1*nperl+i2;

			int i3 = i1+0;
			int j3 = i2+0;

			int i4 = i1+1;
			int j4 = i2+0;

			int i5 = i1-1;
			int j5 = i2+1;

			int i6 = i1+0;
			int j6 = i2+1;

			int i7 = i1+1;
			int j7 = i2+1;

			prdshft(i3,nperl);
			prdshft(j3,nperl);

			prdshft(i4,nperl);
			prdshft(j4,nperl);

			prdshft(i5,nperl);
			prdshft(j5,nperl);
			
			prdshft(i6,nperl);
			prdshft(j6,nperl);
			
			prdshft(i7,nperl);
			prdshft(j7,nperl);		

			cells1[itery] =  b1;
			cells2[itery] =  i3*nperl+j3;

			itery++;

			cells1[itery] =  b1;
			cells2[itery] =  i4*nperl+j4;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i5*nperl+j5;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i6*nperl+j6;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i7*nperl+j7;

			itery++;


		}
	}
	int size4 = nbpairs*sizeof(int);

	int *d_cells1;
	int *d_cells2;

	cudaMalloc((void**)&d_cells1,size4);

	cudaMalloc((void**)&d_cells2,size4);

	cudaMemcpy(d_cells1,cells1,size4,cudaMemcpyHostToDevice);
	cudaMemcpy(d_cells2,cells2,size4,cudaMemcpyHostToDevice);
	//int sibdiv = floor(ll/4.0);
	// print_device_array(d_cells1,nbpairs);
	// print_device_array(d_cells2,nbpairs);



	// matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);
	
	float2 *particles = new float2 [totalN];
	float2 *momenta = new float2 [totalN];
	int *p_indices = new int [totalN];

	for(int i = 0 ; i < totalN ; i++)
	p_indices[i]=i;

	float2 *d_particles;
	float2 *d_momenta;
	int *d_p_indices;

	int *d_bound;
	double *d_boundalong;
	int *d_changestate;

	cudaMalloc((void**)&d_bound,(na+nb)*sizeof(int));
	cudaMalloc((void**)&d_boundalong,(na+nb)*sizeof(double));
	cudaMalloc((void**)&d_changestate,(na+nb)*sizeof(int));

	cudaMemset(d_bound,0,(na+nb)*sizeof(int));
	cudaMemset(d_boundalong,0.,(na+nb)*sizeof(double));
	cudaMemset(d_changestate,0,(na+nb)*sizeof(int));

	double *d_totalforcex;
	double *d_totalforcey;

	cudaMalloc((void**)&d_totalforcex,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey,totalN*sizeof(double));


	matrix<double> state(obj->getdat());


	for(int i = 0  ; i < totalN ; i++) {

	float2 c;
	c.x=state(i,0);
	c.y=state(i,1);

	(particles)[i]=c;

	float2 d;

	d.x = 0.;
	d.y = 0.;

	(momenta)[i]=d;
	}


	int size =  totalN*sizeof(float2);
	int size2 = totalN*sizeof(int);


	cudaMalloc((void**)&d_particles,size);
	cudaMalloc((void**)&d_momenta,size);
	cudaMalloc((void**)&d_p_indices,size2);

	cudaMemcpy(d_particles,particles,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_momenta,momenta,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_p_indices,p_indices,size2,cudaMemcpyHostToDevice);


	// matrix<int> *froyo1 = obj->calculatepairs(boxes,pai,3.5);
	// matrix<int> *froyo2 = obj->calculatepairs(boxes,pbi,3.5);
	// matrix<int> *froyo3 = obj->calculatepairs(boxes,pci,3.5);
	// matrix<int> *froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
	// matrix<int> *froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
	// matrix<int> *froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);
	
	int *d_indices1;
	int *d_indices2;
	double *d_close;


	int tpp;
	


	construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp);



	less_than_condition_AND cond1(SQR(3.5),0,na);
	less_than_condition_AND cond2(SQR(3.5),na,na+nb);
	less_than_condition_AND cond3(SQR(3.5),na+nb,na+nb+nc);
	less_than_condition_NAND cond4(SQR(3.5),0,na,na,na+nb);
	less_than_condition_NAND cond5(SQR(3.5),0,na,na+nb,na+nb+nc);
	less_than_condition_NAND cond6(SQR(3.5),na,na+nb,na+nb,na+nb+nc);


	int th1;
	int *d1_list1,*d1_list2,*d1_list3,*d1_list4;
	pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa

	int th2;
	int *d2_list1,*d2_list2,*d2_list3,*d2_list4;
	pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2);	//fbb

	int th3;
	int *d3_list1,*d3_list2,*d3_list3,*d3_list4;
	pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3);	//fcc

	int th4;
	int *d4_list1,*d4_list2,*d4_list3,*d4_list4;
	pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4);	//fab

	int th5;
	int *d5_list1,*d5_list2,*d5_list3,*d5_list4;
	pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5);	//fac

	int th6;
	int *d6_list1,*d6_list2,*d6_list3,*d6_list4;
	pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc			
	//matrix<double> state(obj->getdat()); //the state of the system

	int th8 = (*bondpairs).getNsafe();
	int th9 = (*bendtriplets).getNsafe();


	int *d8_list1,*d8_list2,*d8_list3,*d8_list4;
	cudaMalloc((void**)&d8_list1,th8*sizeof(int));
	cudaMalloc((void**)&d8_list2,th8*sizeof(int));
	cudaMalloc((void**)&d8_list3,th8*sizeof(int));
	cudaMalloc((void**)&d8_list4,th8*sizeof(int));



	int *h8_list1 = new int [th8];
	int *h8_list2 = new int [th8];

	for(int i = 0 ; i < th8 ; i++ ) {
		h8_list1[i] = (*bondpairs)(i,0);
		h8_list2[i] = (*bondpairs)(i,1);
	}

	cudaMemcpy(d8_list1,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list2,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list3,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list4,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);





	int *d9_list1,*d9_list2,*d9_list3,*d9_list4,*d9_list5,*d9_list6;
	cudaMalloc((void**)&d9_list1,th9*sizeof(int));
	cudaMalloc((void**)&d9_list2,th9*sizeof(int));
	cudaMalloc((void**)&d9_list3,th9*sizeof(int));
	cudaMalloc((void**)&d9_list4,th9*sizeof(int));
	cudaMalloc((void**)&d9_list5,th9*sizeof(int));
	cudaMalloc((void**)&d9_list6,th9*sizeof(int));

	int *h9_list1 = new int [th9];
	int *h9_list2 = new int [th9];
	int *h9_list3 = new int [th9];

	for(int i = 0 ; i < th9 ; i++ ) {
		h9_list1[i] = (*bendtriplets)(i,0);
		h9_list2[i] = (*bendtriplets)(i,1);
		h9_list3[i] = (*bendtriplets)(i,2);
	}	

	cudaMemcpy(d9_list1,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list2,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list3,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list4,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list5,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list6,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);


	int i;

	double cons1;
	double cons2;
	double cons3;
	double cons4;

	//(mom)->operator()(i,i1) = c5*c2*((mom)->operator()(i,i1)) + (c5*(c3)+q)*F(i,i1) + (c5*(c4)+r)*R(i,i1);
	cons1 = (*obj).getc5()*(*obj).getc2();
	cons2 = (*obj).getc5()*(*obj).getc3()+(*obj).getq();
	cons3 = (*obj).getc5()*(*obj).getc4()+(*obj).getr();
	cons4 = (*obj).getc1();

	double d_dt = (*obj).getdt();
	double d_m = (*obj).getm();
	double d_kT = (*obj).getkT();


	// vector<matrix<double> > savef1forces;
	// vector<matrix<double> > savef2forces;
	// vector<matrix<double> > savef3forces;
	// vector<matrix<double> > savef4forces;
	// vector<matrix<double> > savef5forces;
	// vector<matrix<double> > savef6forces;
	// vector<matrix<double> > savef7forces;

	// vector<matrix<double> > saveftemp1forces;
	// vector<matrix<double> > saveftemp2forces;
	// vector<matrix<double> > saveftemp3forces;

	// vector<matrix<double> > savepositions;
	// vector<vector1<int> > savebound;
	// vector<vector1<double> > saveboundalongs;

	for(i = 0 ; i < runtime ; i++) {
		//cout << i << endl;

		cout << i << endl;
	
		//cout << (*obj).avmom() << endl;
	if(i%25==0) {
		//delete froyo1,froyo2,froyo3,froyo4,froyo5,froyo6;
		// cout << "updated after: " << i << endl;
		// state = obj->getdat();


		cudaFree(d1_list1);cudaFree(d1_list2);cudaFree(d1_list3);cudaFree(d1_list4);
		

		cudaFree(d2_list1);cudaFree(d2_list2);cudaFree(d2_list3);cudaFree(d2_list4);
		

		cudaFree(d3_list1);cudaFree(d3_list2);cudaFree(d3_list3);cudaFree(d3_list4);
		

		cudaFree(d4_list1);cudaFree(d4_list2);cudaFree(d4_list3);cudaFree(d4_list4);
		

		cudaFree(d5_list1);cudaFree(d5_list2);cudaFree(d5_list3);cudaFree(d5_list4);
	
		cudaFree(d6_list1);cudaFree(d6_list2);cudaFree(d6_list3);cudaFree(d6_list4);

		cudaFree(d_indices1);

		cudaFree(d_indices2);

		cudaFree(d_close);

		this->resetindices(d_p_indices,totalN);


		construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp,false);

		// cout << tpp << endl;
		// cout << "pair" << endl;

		pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa
		pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2); //fbb
		pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3); //fcc
		pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4); //fab
		pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5); //fac
		pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc

	}




	double *d1_forces1x;
	double *d1_forces2x;
	double *d1_forces1y;
	double *d1_forces2y;
	cudaMalloc((void**)&d1_forces1x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces1y,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2y,th1*sizeof(double));
	calculateforces2D(d1_list1,d1_list2,d_particles, d1_forces1x,d1_forces1y,d1_forces2x,d1_forces2y, faa_gpu ,th1, l,true);

	// arracychck(d1_forces1x,th1);
	// arracychck(d1_forces1y,th1);
	// cout << "force1" << endl;

	double *d2_forces1x;
	double *d2_forces2x;
	double *d2_forces1y;
	double *d2_forces2y;
	cudaMalloc((void**)&d2_forces1x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces1y,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2y,th2*sizeof(double));
	calculateforces2D(d2_list1,d2_list2,d_particles, d2_forces1x,d2_forces1y,d2_forces2x,d2_forces2y, fbb_gpu ,th2, l,true);	



	// arracychck(d2_forces1x,th2);
	// arracychck(d2_forces1y,th2);

	// cout << "force2" << endl;
	double *d3_forces1x;
	double *d3_forces2x;
	double *d3_forces1y;
	double *d3_forces2y;
	cudaMalloc((void**)&d3_forces1x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces1y,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2y,th3*sizeof(double));
	calculateforces2D(d3_list1,d3_list2,d_particles, d3_forces1x,d3_forces1y,d3_forces2x,d3_forces2y, fcc_gpu ,th3, l,true);	

	// arracychck(d3_forces1x,th3);
	// arracychck(d3_forces1y,th3);

	// cout << "force3" << endl;
	double *d4_forces1x;
	double *d4_forces2x;
	double *d4_forces1y;
	double *d4_forces2y;
	cudaMalloc((void**)&d4_forces1x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces1y,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2y,th4*sizeof(double));
	calculateforces2D(d4_list1,d4_list2,d_particles, d4_forces1x,d4_forces1y,d4_forces2x,d4_forces2y, fab_gpu ,th4, l,true);
	
	// arracychck(d4_forces1x,th4);
	// arracychck(d4_forces1y,th4);
	// cout << "force4" << endl;

	double *d5_forces1x;
	double *d5_forces2x;
	double *d5_forces1y;
	double *d5_forces2y;
	cudaMalloc((void**)&d5_forces1x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces1y,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2y,th5*sizeof(double));
	calculateforces2D(d5_list1,d5_list2,d_particles, d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y, fbc_gpu ,th5, l,true);
	// arracychck(d5_forces1x,th5);
	// arracychck(d5_forces1y,th5);
	// cout << "force5" << endl;
	// cout << "d5" << endl;
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y,d_particles,th5,totalN);
	// pausel();	

	double *d6_forces1x;
	double *d6_forces2x;
	double *d6_forces1y;
	double *d6_forces2y;
	cudaMalloc((void**)&d6_forces1x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces1y,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2y,th6*sizeof(double));
	calculateforces2D(d6_list1,d6_list2,d_particles, d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y, fac_gpu ,th6, l,true);	
	// arracychck(d6_forces1x,th6);
	// arracychck(d6_forces1y,th6);
	// cout << "force6" << endl;
	// cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();

	// matrix<double> F6((*obj).calculateforces(*bondpairs,*bindm));
	callCalculateUnbindingsGPU(d_particles,d_bound,d_boundalong,d_changestate);



	// cout << "unbindings calculated" << endl;


	callCalculateBindingsGPU(d5_list1,d5_list2,d6_list1,d6_list2,d_particles, d_bound, d_boundalong,d_changestate,th5 ,th6 );


	// cout << "bindings calculated" << endl;


	int *d7_list1,*d7_list2,*d7_list3;
	double *d7_forces1x;
	double *d7_forces1y;
	double *d7_forces2x;
	double *d7_forces2y;
	double *d7_forces3x;
	double *d7_forces3y;
	int th7;
	BindingForcesGPU(d_particles, d_bound, d_boundalong, d7_list1,d7_list2,d7_list3, d7_forces1x, d7_forces1y, d7_forces2x,d7_forces2y,d7_forces3x, d7_forces3y, bindp_gpu, th7);

	// arracychck(d7_forces1x,th7);
	// arracychck(d7_forces1y,th7);
	// cout << "binding forces calculated" << endl;


	double *d8_forces1x;
	double *d8_forces2x;
	double *d8_forces1y;
	double *d8_forces2y;
	cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));
	calculateforces2D(d8_list1,d8_list2,d_particles, d8_forces1x,d8_forces1y,d8_forces2x,d8_forces2y, bindm_gpu ,th8, l,true);
	
	// arracychck(d8_forces1x,th8);
	// arracychck(d8_forces1y,th8);

	// cout << "force 8" << endl;
	resetchangestate(d_changestate);
	// matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets,*bendp));

//	int th9;


	double *d9_forces1x;
	double *d9_forces2x;
	double *d9_forces1y;
	double *d9_forces2y;
	double *d9_forces3x;
	double *d9_forces3y;
	// cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));	
	cudaMalloc((void**)&d9_forces1x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces1y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3y,th9*sizeof(double));
	BendingForcesGPU(d_particles, d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces1y,d9_forces2x,d9_forces2y, d9_forces3x, d9_forces3y,bendp_gpu,th9);

	// print_device_array_weave(d9_forces1x,d9_forces1y,th9);
	// print_device_array_weave(d9_forces2x,d9_forces2y,th9);
	// print_device_array_weave(d9_forces3x,d9_forces3y,th9);
	// cout << "force 9" << endl;	
	// arracychck(d9_forces1x,th9);
	// arracychck(d9_forces1y,th9);

	// print_device_float2(d_particles,totalN);
	// print_device_array(d9_list1,th9);
	// print_device_array(d9_list2,th9);
	// print_device_array(d9_list3,th9);

	// cout << "force9" << endl;

	//	matrix<double> F = ftemp1+ftemp2+ftemp3+F1+F2+F3+F4+F5+F6+F7;//+F4+F5;

		// matrix<double> R(totalN,dimension);
		// for(int i1 = 0 ; i1 < totalN ; i1++) {
		// 	for(int j = 0 ; j < dimension ; j++) {
		// 		R(i1,j) = (3.464101615 * ((double) rand() / (RAND_MAX)) - 1.732050808);
		// 	}
		// }

	int *d10_list1;
	double *d10_forces1x;
	double *d10_forces1y;	
	int th10;

	PositionForcesDueToAnglesGPU(d_particles, d_bound, d_boundalong, d10_list1, d10_forces1x, d10_forces1y,th10);
	

	// arracychck(d10_forces1x,th10);
	// arracychck(d10_forces1y,th10);
	// cout << "force10" << endl;
	// cout << "all forces calculated" << endl;





	resetforce(d_totalforcex);

	resetforce(d_totalforcey);

	//cout << "reset" << endl;

	// print_device_array(d_totalforcex,totalN);


	ReduceForces(d1_list1,d1_list2,d1_list3,d1_list4,d1_forces1x,d1_forces2x,d1_forces1y,d1_forces2y,d_totalforcex,d_totalforcey,th1);
	// cout << "d1" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d2_list1,d2_list2,d2_list3,d2_list4,d2_forces1x,d2_forces2x,d2_forces1y,d2_forces2y,d_totalforcex,d_totalforcey,th2);
	/// cout << "d2" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces(d3_list1,d3_list2,d3_list3,d3_list4,d3_forces1x,d3_forces2x,d3_forces1y,d3_forces2y,d_totalforcex,d_totalforcey,th3);
	// cout << "d3" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d4_list1,d4_list2,d4_list3,d4_list4,d4_forces1x,d4_forces2x,d4_forces1y,d4_forces2y,d_totalforcex,d_totalforcey,th4);
	// cout << "d4" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	
	ReduceForces(d5_list1,d5_list2,d5_list3,d5_list4,d5_forces1x,d5_forces2x,d5_forces1y,d5_forces2y,d_totalforcex,d_totalforcey,th5);
	// cout << "d5" << endl;	
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d_particles,th5,totalN);
	// pausel();	
	ReduceForces(d6_list1,d6_list2,d6_list3,d6_list4,d6_forces1x,d6_forces2x,d6_forces1y,d6_forces2y,d_totalforcex,d_totalforcey,th6);	
	// cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();
	ReduceForces3(d7_list1,d7_list2,d7_list3,d7_forces1x,d7_forces2x,d7_forces3x,d7_forces1y,d7_forces2y,d7_forces3y,d_totalforcex,d_totalforcey,th7);	
	// cout << "d7" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces(d8_list1,d8_list2,d8_forces1x,d8_forces2x,d8_forces1y,d8_forces2y,d_totalforcex,d_totalforcey,th8);
	// cout << "d8" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	ReduceForces3(d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces2x,d9_forces3x,d9_forces1y,d9_forces2y,d9_forces3y,d_totalforcex,d_totalforcey,th9);		
	
	// cout << "d9" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double ff = (v0_a+v0_b)/2.;
	ReduceForcesAndNormalize(d10_list1,d10_forces1x,d10_forces1y,d_totalforcex,d_totalforcey, max_s, ff, th10);

	cout << "reduction" << endl;
	// cout << "d10" << endl;
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	if(i>0&&i%every==0) { 
		// for(int j = 0 ; j < na+nb ; j++) {
		// if(bound[j]>0){
		// cout << "printed" << endl;
		// cout << j << endl;
		// cout << bound[j] << endl;
		// cout << bound_along[j] << endl;
		// cout << F5(j,'r') << endl;
		// cout << F4(j,'r') << endl;
		// cout << obj->getcoordinate(j,0) << " " << obj->getcoordinate(j,1) << endl;
		// cout << obj->getcoordinate(200,0) <<  " " << obj->getcoordinate(200,1) << endl;
		// cout << obj->getcoordinate(200+L,0) <<  " " << obj->getcoordinate(200+L,1) << endl;
		// }
		// }
		// cout << F6 << endl;
		// cout << F7 << endl;

		// cout << ftemp2 << endl;
		// cout << ftemp3 << endl;

		stringstream ss2;
		// ss2 <<i/every;
		// string pairlist = "list";

		 stringstream kts;
		 kts << (*obj).getkT();

		 // stringstream epi;
		 // epi << eps;

		 // stringstream epieq;
		 // epieq << eqeps;

		 stringstream len;
		 len << l;

		 string extension =  "_kT="+kts.str()+"_l="+len.str()+".csv";

		stringstream ss;
		ss <<(i/every);
		string filename = "x";
		filename += ss.str();
		filename += extension;

		string momname = "bind";
		momname += ss.str();
		momname += extension;

		string baname = "bind_along";
		baname += ss.str();
		baname += extension;



		ofstream myfile;
		myfile.open(filename.c_str());
		//myfile <<= (*obj).getdat();
		file_print_device_float2(d_particles,totalN,myfile);
		myfile.close();



		}	

	double *d_R1;
	double *d_R2;

	// cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));	
	cudaMalloc((void**)&d_R1,totalN*sizeof(double));
	cudaMalloc((void**)&d_R2,totalN*sizeof(double));

	setstaterandom(d_R1,1.732050808,totalN);
 	setstaterandom(d_R2,1.732050808,totalN);	

 	//advmom2D(d_momenta, d_totalforcex, d_totalforcey, d_R1, d_R2, cons1,cons2,cons3,totalN);
 	advmom2D_particledependence(d_momenta,d_totalforcex,d_totalforcey,d_R1,d_R2,func, d_dt, d_kT, d_m, totalN);
 	advpos2D(d_particles, d_momenta, cons4, totalN);


 	applypbc2D(d_particles,d_momenta,l,is_periodic,totalN);

 	cout << "updated" << endl;

	cudaFree(d1_forces1x);
	cudaFree(d1_forces2x);
	cudaFree(d1_forces1y);
	cudaFree(d1_forces2y);
	cudaFree(d2_forces1x);
	cudaFree(d2_forces2x);
	cudaFree(d2_forces1y);
	cudaFree(d2_forces2y);
	cudaFree(d3_forces1x);
	cudaFree(d3_forces2x);
	cudaFree(d3_forces1y);
	cudaFree(d3_forces2y);
	cudaFree(d4_forces1x);
	cudaFree(d4_forces2x);
	cudaFree(d4_forces1y);
	cudaFree(d4_forces2y);
	cudaFree(d5_forces1x);
	cudaFree(d5_forces2x);
	cudaFree(d5_forces1y);
	cudaFree(d5_forces2y);
	cudaFree(d6_forces1x);
	cudaFree(d6_forces2x);
	cudaFree(d6_forces1y);
	cudaFree(d6_forces2y);
 	cudaFree(d7_list1);
 	cudaFree(d7_list2);
 	cudaFree(d7_list3);
	cudaFree(d7_forces1x);
	cudaFree(d7_forces2x);
	cudaFree(d7_forces3x);
	cudaFree(d7_forces1y);
	cudaFree(d7_forces2y);
	cudaFree(d7_forces3y);
	cudaFree(d8_forces1x);
	cudaFree(d8_forces2x);
	cudaFree(d8_forces1y);
	cudaFree(d8_forces2y);
	cudaFree(d9_forces1x);
	cudaFree(d9_forces2x);
	cudaFree(d9_forces3x);
	cudaFree(d9_forces1y);
	cudaFree(d9_forces2y);
	cudaFree(d9_forces3y);
	cudaFree(d10_list1);
	cudaFree(d10_forces1x);
	cudaFree(d10_forces1y);

	cout << "freed" << endl;


	}
}
/*
template <typename Fun>
void Microtubule::runGPUcheck(int runtime, int every, Fun func)
{
//	pausel();
	int ccc;
	int totalN = obj->getN();

	//num is the number of boxes per length

	int ncells = num*num;

	WCApotentialGPU faa_gpu(2.,1.,2.);
	WCApotentialGPU fab_gpu(1.,1.,0.);
	WCApotentialGPU fac_gpu(1.,1.,0.);
	WCApotentialGPU fbb_gpu(2.,1.,2.);
	WCApotentialGPU fbc_gpu(1.,1.,0.);
	WCApotentialGPU fcc_gpu(1.,1.,0.);
	HarmonicPotentialGPU bindp_gpu(100.,0.); 
	FENEPotentialGPU bindm_gpu(50.,1.5); 
	BendingPotentialGPU bendp_gpu(100.,0.);

//we now have the count of each cell list

	matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);


	matrix<int> *froyo1 = obj->calculatepairs(boxes,pai,3.5);
	matrix<int> *froyo2 = obj->calculatepairs(boxes,pbi,3.5);
	matrix<int> *froyo3 = obj->calculatepairs(boxes,pci,3.5);
	matrix<int> *froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
	matrix<int> *froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
	matrix<int> *froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);


	int nbpairs = 5*ncells;
	int nperl = num;

	int *cells1 = new int [nbpairs];
	int *cells2 = new int [nbpairs];


	int itery = 0;
	for(int i1 = 0 ; i1 < num ; i1++) {
		for(int i2 = 0 ; i2 < num ; i2++ ) {


			int b1 =  i1*nperl+i2;

			int i3 = i1+0;
			int j3 = i2+0;

			int i4 = i1+1;
			int j4 = i2+0;

			int i5 = i1-1;
			int j5 = i2+1;

			int i6 = i1+0;
			int j6 = i2+1;

			int i7 = i1+1;
			int j7 = i2+1;

			prdshft(i3,nperl);
			prdshft(j3,nperl);

			prdshft(i4,nperl);
			prdshft(j4,nperl);

			prdshft(i5,nperl);
			prdshft(j5,nperl);
			
			prdshft(i6,nperl);
			prdshft(j6,nperl);
			
			prdshft(i7,nperl);
			prdshft(j7,nperl);		

			cells1[itery] =  b1;
			cells2[itery] =  i3*nperl+j3;

			itery++;

			cells1[itery] =  b1;
			cells2[itery] =  i4*nperl+j4;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i5*nperl+j5;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i6*nperl+j6;
			
			itery++;
			
			cells1[itery] =  b1;
			cells2[itery] =  i7*nperl+j7;

			itery++;


		}
	}
	int size4 = nbpairs*sizeof(int);

	int *d_cells1;
	int *d_cells2;

	cudaMalloc((void**)&d_cells1,size4);

	cudaMalloc((void**)&d_cells2,size4);

	cudaMemcpy(d_cells1,cells1,size4,cudaMemcpyHostToDevice);
	cudaMemcpy(d_cells2,cells2,size4,cudaMemcpyHostToDevice);
	//int sibdiv = floor(ll/4.0);
	// print_device_array(d_cells1,nbpairs);
	// print_device_array(d_cells2,nbpairs);



	// matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);
	
	float2 *particles = new float2 [totalN];
	float2 *momenta = new float2 [totalN];
	int *p_indices = new int [totalN];

	for(int i = 0 ; i < totalN ; i++)
	p_indices[i]=i;

	float2 *d_particles;
	float2 *d_momenta;
	int *d_p_indices;

	int *d_bound;
	double *d_boundalong;
	int *d_changestate;

	cudaMalloc((void**)&d_bound,(na+nb)*sizeof(int));
	cudaMalloc((void**)&d_boundalong,(na+nb)*sizeof(double));
	cudaMalloc((void**)&d_changestate,(na+nb)*sizeof(int));

	cudaMemset(d_bound,0,(na+nb)*sizeof(int));
	cudaMemset(d_boundalong,0.,(na+nb)*sizeof(double));
	cudaMemset(d_changestate,0,(na+nb)*sizeof(int));

	double *d_totalforcex;
	double *d_totalforcey;

	cudaMalloc((void**)&d_totalforcex,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey,totalN*sizeof(double));


	matrix<double> state(obj->getdat());


	for(int i = 0  ; i < totalN ; i++) {

	float2 c;
	c.x=state(i,0);
	c.y=state(i,1);

	(particles)[i]=c;

	float2 d;

	d.x = 0.;
	d.y = 0.;

	(momenta)[i]=d;
	}


	int size =  totalN*sizeof(float2);
	int size2 = totalN*sizeof(int);


	cudaMalloc((void**)&d_particles,size);
	cudaMalloc((void**)&d_momenta,size);
	cudaMalloc((void**)&d_p_indices,size2);

	cudaMemcpy(d_particles,particles,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_momenta,momenta,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_p_indices,p_indices,size2,cudaMemcpyHostToDevice);


	// matrix<int> *froyo1 = obj->calculatepairs(boxes,pai,3.5);
	// matrix<int> *froyo2 = obj->calculatepairs(boxes,pbi,3.5);
	// matrix<int> *froyo3 = obj->calculatepairs(boxes,pci,3.5);
	// matrix<int> *froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
	// matrix<int> *froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
	// matrix<int> *froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);
	
	int *d_indices1;
	int *d_indices2;
	double *d_close;


	int tpp;
	


	construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp);



	less_than_condition_AND cond1(SQR(3.5),0,na);
	less_than_condition_AND cond2(SQR(3.5),na,na+nb);
	less_than_condition_AND cond3(SQR(3.5),na+nb,na+nb+nc);
	less_than_condition_NAND cond4(SQR(3.5),0,na,na,na+nb);
	less_than_condition_NAND cond5(SQR(3.5),0,na,na+nb,na+nb+nc);
	less_than_condition_NAND cond6(SQR(3.5),na,na+nb,na+nb,na+nb+nc);


	int th1;
	int *d1_list1,*d1_list2,*d1_list3,*d1_list4;
	pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa

	int th2;
	int *d2_list1,*d2_list2,*d2_list3,*d2_list4;
	pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2);	//fbb

	int th3;
	int *d3_list1,*d3_list2,*d3_list3,*d3_list4;
	pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3);	//fcc

	int th4;
	int *d4_list1,*d4_list2,*d4_list3,*d4_list4;
	pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4);	//fab

	int th5;
	int *d5_list1,*d5_list2,*d5_list3,*d5_list4;
	pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5);	//fac

	int th6;
	int *d6_list1,*d6_list2,*d6_list3,*d6_list4;
	pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc			
	//matrix<double> state(obj->getdat()); //the state of the system

	int th8 = (*bondpairs).getNsafe();
	int th9 = (*bendtriplets).getNsafe();


	int *d8_list1,*d8_list2,*d8_list3,*d8_list4;
	cudaMalloc((void**)&d8_list1,th8*sizeof(int));
	cudaMalloc((void**)&d8_list2,th8*sizeof(int));
	cudaMalloc((void**)&d8_list3,th8*sizeof(int));
	cudaMalloc((void**)&d8_list4,th8*sizeof(int));



	int *h8_list1 = new int [th8];
	int *h8_list2 = new int [th8];

	for(int i = 0 ; i < th8 ; i++ ) {
		h8_list1[i] = (*bondpairs)(i,0);
		h8_list2[i] = (*bondpairs)(i,1);
	}

	cudaMemcpy(d8_list1,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list2,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list3,h8_list1,th8*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d8_list4,h8_list2,th8*sizeof(int),cudaMemcpyHostToDevice);





	int *d9_list1,*d9_list2,*d9_list3,*d9_list4,*d9_list5,*d9_list6;
	cudaMalloc((void**)&d9_list1,th9*sizeof(int));
	cudaMalloc((void**)&d9_list2,th9*sizeof(int));
	cudaMalloc((void**)&d9_list3,th9*sizeof(int));
	cudaMalloc((void**)&d9_list4,th9*sizeof(int));
	cudaMalloc((void**)&d9_list5,th9*sizeof(int));
	cudaMalloc((void**)&d9_list6,th9*sizeof(int));

	int *h9_list1 = new int [th9];
	int *h9_list2 = new int [th9];
	int *h9_list3 = new int [th9];

	for(int i = 0 ; i < th9 ; i++ ) {
		h9_list1[i] = (*bendtriplets)(i,0);
		h9_list2[i] = (*bendtriplets)(i,1);
		h9_list3[i] = (*bendtriplets)(i,2);
	}	

	cudaMemcpy(d9_list1,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list2,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list3,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list4,h9_list1,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list5,h9_list2,th9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d9_list6,h9_list3,th9*sizeof(int),cudaMemcpyHostToDevice);


	int i;

	double cons1;
	double cons2;
	double cons3;
	double cons4;

	//(mom)->operator()(i,i1) = c5*c2*((mom)->operator()(i,i1)) + (c5*(c3)+q)*F(i,i1) + (c5*(c4)+r)*R(i,i1);
	cons1 = (*obj).getc5()*(*obj).getc2();
	cons2 = (*obj).getc5()*(*obj).getc3()+(*obj).getq();
	cons3 = (*obj).getc5()*(*obj).getc4()+(*obj).getr();
	cons4 = (*obj).getc1();


	// vector<matrix<double> > savef1forces;
	// vector<matrix<double> > savef2forces;
	// vector<matrix<double> > savef3forces;
	// vector<matrix<double> > savef4forces;
	// vector<matrix<double> > savef5forces;
	// vector<matrix<double> > savef6forces;
	// vector<matrix<double> > savef7forces;

	// vector<matrix<double> > saveftemp1forces;
	// vector<matrix<double> > saveftemp2forces;
	// vector<matrix<double> > saveftemp3forces;

	// vector<matrix<double> > savepositions;
	// vector<vector1<int> > savebound;
	// vector<vector1<double> > saveboundalongs;

	for(i = 0 ; i < runtime ; i++) {
		//cout << i << endl;

		cout << i << endl;
		//pausel();
	
		//cout << (*obj).avmom() << endl;
	if(i%25==0) {
		//delete froyo1,froyo2,froyo3,froyo4,froyo5,froyo6;
		// cout << "updated after: " << i << endl;
		// state = obj->getdat();

		delete froyo1,froyo2,froyo3,froyo4,froyo5,froyo6;
		// cout << "updated after: " << i << endl;
		// state = obj->getdat();
		froyo1 = obj->calculatepairs(boxes,pai,3.5);
		froyo2 = obj->calculatepairs(boxes,pbi,3.5);
		froyo3 = obj->calculatepairs(boxes,pci,3.5);
		froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
		froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
		froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);


		cudaFree(d1_list1);cudaFree(d1_list2);cudaFree(d1_list3);cudaFree(d1_list4);
		

		cudaFree(d2_list1);cudaFree(d2_list2);cudaFree(d2_list3);cudaFree(d2_list4);
		

		cudaFree(d3_list1);cudaFree(d3_list2);cudaFree(d3_list3);cudaFree(d3_list4);
		

		cudaFree(d4_list1);cudaFree(d4_list2);cudaFree(d4_list3);cudaFree(d4_list4);
		

		cudaFree(d5_list1);cudaFree(d5_list2);cudaFree(d5_list3);cudaFree(d5_list4);
	
		cudaFree(d6_list1);cudaFree(d6_list2);cudaFree(d6_list3);cudaFree(d6_list4);

		cudaFree(d_indices1);

		cudaFree(d_indices2);

		cudaFree(d_close);



		this->resetindices(d_p_indices,totalN);



		construct_possible_pair_list(d_particles,d_p_indices,totalN,l,d_cells1,d_cells2,num,is_periodic,d_indices1,d_indices2,d_close,tpp,false);


		cout << "pair" << endl;

		pairlist(d_indices1,d_indices2,d_close,cond1,d1_list1,d1_list2,d1_list3,d1_list4,tpp, th1); //faa
		pairlist(d_indices1,d_indices2,d_close,cond2,d2_list1,d2_list2,d2_list3,d2_list4,tpp, th2); //fbb
		pairlist(d_indices1,d_indices2,d_close,cond3,d3_list1,d3_list2,d3_list3,d3_list4,tpp, th3); //fcc
		pairlist(d_indices1,d_indices2,d_close,cond4,d4_list1,d4_list2,d4_list3,d4_list4,tpp, th4); //fab
		pairlist(d_indices1,d_indices2,d_close,cond5,d5_list1,d5_list2,d5_list3,d5_list4,tpp, th5); //fac
		pairlist(d_indices1,d_indices2,d_close,cond6,d6_list1,d6_list2,d6_list3,d6_list4,tpp, th6); //fbc


		
		// froyo1 = obj->calculatepairs(boxes,pai,3.5);
		// froyo2 = obj->calculatepairs(boxes,pbi,3.5);
		// froyo3 = obj->calculatepairs(boxes,pci,3.5);
		// froyo4 = obj->calculatepairs(boxes,pai,pbi,3.5);
		// froyo5 = obj->calculatepairs(boxes,pai,pci,3.5);
		// froyo6 = obj->calculatepairs(boxes,pbi,pci,3.5);

	}





	//cout << "pairs" << endl;
	// cout << "pairings" << endl;


	// matrix<double> ftemp2(totalN,dimension),ftemp3(totalN,dimension);
	// //matrix<double> angforces1(nc,dimension-1),angforces2(nc,dimension-1),angforces3(nc,dimension-1);
	// // cout << "matrices initialized" << endl;

	// matrix<double> F1((*obj).calculateforces(*froyo1,*faa)); //calculate the forces using the pairs as an input

	// matrix<double> F2((*obj).calculateforces(*froyo2,*fbb)); //calculate the forces using the pairs as an input

	// matrix<double> ftemp1((*obj).calculateforces(*froyo3,*fcc)); //calculate the forces using the pairs as an input
	
	// matrix<double> F3((*obj).calculateforces(*froyo4,*fab)); //calculate the forces using the pairs as an input

	// this->ForcesDueToPositionPL(*froyo5,ftemp2); //calculate the forces using the pairs as an input

	// this->ForcesDueToPositionPL(*froyo6,ftemp3); //calculate the forces using the pairs as an input

	// this->CalculateBindings(*froyo5,*froyo6);

	// matrix<double> F4 = this->BindingForces();

	// matrix<double> F5 = this->PositionForcesDueToAngles();

	//print_device_float2(d_particles,totalN);

	double *d1_forces1x;
	double *d1_forces2x;
	double *d1_forces1y;
	double *d1_forces2y;
	cudaMalloc((void**)&d1_forces1x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces1y,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2x,th1*sizeof(double));
	cudaMalloc((void**)&d1_forces2y,th1*sizeof(double));
	calculateforces2D(d1_list1,d1_list2,d_particles, d1_forces1x,d1_forces1y,d1_forces2x,d1_forces2y, faa_gpu ,th1, l,true);

	cout << "force1" << endl;
	arracychck(d1_forces1x,th1);
	arracychck(d1_forces1y,th1); 

	double *d2_forces1x;
	double *d2_forces2x;
	double *d2_forces1y;
	double *d2_forces2y;
	cudaMalloc((void**)&d2_forces1x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces1y,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2x,th2*sizeof(double));
	cudaMalloc((void**)&d2_forces2y,th2*sizeof(double));
	calculateforces2D(d2_list1,d2_list2,d_particles, d2_forces1x,d2_forces1y,d2_forces2x,d2_forces2y, fbb_gpu ,th2, l,true);	

	cout << "force2" << endl;
	arracychck(d1_forces2x,th2);
	arracychck(d1_forces2y,th2); 

	double *d3_forces1x;
	double *d3_forces2x;
	double *d3_forces1y;
	double *d3_forces2y;
	cudaMalloc((void**)&d3_forces1x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces1y,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2x,th3*sizeof(double));
	cudaMalloc((void**)&d3_forces2y,th3*sizeof(double));
	calculateforces2D(d3_list1,d3_list2,d_particles, d3_forces1x,d3_forces1y,d3_forces2x,d3_forces2y, fcc_gpu ,th3, l,true);	


	cout << "force3" << endl;
	arracychck(d3_forces1x,th3);
	arracychck(d3_forces1y,th3); 

	double *d4_forces1x;
	double *d4_forces2x;
	double *d4_forces1y;
	double *d4_forces2y;
	cudaMalloc((void**)&d4_forces1x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces1y,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2x,th4*sizeof(double));
	cudaMalloc((void**)&d4_forces2y,th4*sizeof(double));
	calculateforces2D(d4_list1,d4_list2,d_particles, d4_forces1x,d4_forces1y,d4_forces2x,d4_forces2y, fab_gpu ,th4, l,true);

	cout << "force4" << endl;
	arracychck(d4_forces1x,th4);
	arracychck(d4_forces1y,th4); 

	double *d5_forces1x;
	double *d5_forces2x;
	double *d5_forces1y;
	double *d5_forces2y;
	cudaMalloc((void**)&d5_forces1x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces1y,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2x,th5*sizeof(double));
	cudaMalloc((void**)&d5_forces2y,th5*sizeof(double));
	calculateforces2D(d5_list1,d5_list2,d_particles, d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y, fbc_gpu ,th5, l,true);

	cout << "force5" << endl;
	arracychck(d5_forces1x,th5);
	arracychck(d5_forces1y,th5); 	
	// cout << "d5" << endl;
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d5_forces2x,d5_forces2y,d_particles,th5,totalN);
	// pausel();	

	double *d6_forces1x;
	double *d6_forces2x;
	double *d6_forces1y;
	double *d6_forces2y;
	cudaMalloc((void**)&d6_forces1x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces1y,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2x,th6*sizeof(double));
	cudaMalloc((void**)&d6_forces2y,th6*sizeof(double));
	calculateforces2D(d6_list1,d6_list2,d_particles, d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y, fac_gpu ,th6, l,true);	

	cout << "force6" << endl;
	arracychck(d6_forces1x,th6);
	arracychck(d6_forces1y,th6); 	
	// cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d6_forces2x,d6_forces2y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();

	// matrix<double> F6((*obj).calculateforces(*bondpairs,*bindm));
	callCalculateUnbindingsGPU(d_particles,d_bound,d_boundalong,d_changestate);



	cout << "unbindings calculated" << endl;




	callCalculateBindingsGPU(d5_list1,d5_list2,d6_list1,d6_list2,d_particles, d_bound, d_boundalong,d_changestate,th5 ,th6 );


	cout << "bindings calculated" << endl;



	int *d7_list1,*d7_list2,*d7_list3;
	double *d7_forces1x;
	double *d7_forces1y;
	double *d7_forces2x;
	double *d7_forces2y;
	double *d7_forces3x;
	double *d7_forces3y;
	int th7;
	BindingForcesGPU(d_particles, d_bound, d_boundalong, d7_list1,d7_list2,d7_list3, d7_forces1x, d7_forces1y, d7_forces2x,d7_forces2y,d7_forces3x, d7_forces3y, bindp_gpu, th7);


	cout << "binding forces calculated" << endl;
	arracychck(d7_forces1x,th7);
	arracychck(d7_forces1y,th7); 

	double *d8_forces1x;
	double *d8_forces2x;
	double *d8_forces1y;
	double *d8_forces2y;
	cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));
	calculateforces2D(d8_list1,d8_list2,d_particles, d8_forces1x,d8_forces1y,d8_forces2x,d8_forces2y, bindm_gpu ,th8, l,true);

	cout << "forces8" << endl;
	arracychck(d8_forces1x,th8);
	arracychck(d8_forces1y,th8); 

	resetchangestate(d_changestate);
	// matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets,*bendp));

//	int th9;


	double *d9_forces1x;
	double *d9_forces2x;
	double *d9_forces1y;
	double *d9_forces2y;
	double *d9_forces3x;
	double *d9_forces3y;
	cudaMalloc((void**)&d9_forces1x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces1y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces2y,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3x,th9*sizeof(double));
	cudaMalloc((void**)&d9_forces3y,th9*sizeof(double));
	BendingForcesGPU(d_particles, d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces1y,d9_forces2x,d9_forces2y, d9_forces3x, d9_forces3y,bendp_gpu,th9);

	cout << "bending forces" << endl;
	arracychck(d9_forces1x,th9);
	arracychck(d9_forces1y,th9); 
	//	matrix<double> F = ftemp1+ftemp2+ftemp3+F1+F2+F3+F4+F5+F6+F7;//+F4+F5;

		// matrix<double> R(totalN,dimension);
		// for(int i1 = 0 ; i1 < totalN ; i1++) {
		// 	for(int j = 0 ; j < dimension ; j++) {
		// 		R(i1,j) = (3.464101615 * ((double) rand() / (RAND_MAX)) - 1.732050808);
		// 	}
		// }

	int *d10_list1;
	double *d10_forces1x;
	double *d10_forces1y;	
	int th10;

	PositionForcesDueToAnglesGPU(d_particles, d_bound, d_boundalong, d10_list1, d10_forces1x, d10_forces1y,th10);
	

	cout << "pos forces" << endl;
	arracychck(d10_forces1x,th9);
	arracychck(d10_forces1y,th9); 

	cout << "all forces calculated" << endl;




	double *d_totalforcex1;
	double *d_totalforcey1;

	cudaMalloc((void**)&d_totalforcex1,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey1,totalN*sizeof(double));

	resetforce(d_totalforcex1);
	resetforce(d_totalforcey1);

	cout << "reset" << endl;

	// print_device_array(d_totalforcex,totalN);


	//ReduceForces(d1_list1,d1_list2,d1_list3,d1_list4,d1_forces1x,d1_forces2x,d1_forces1y,d1_forces2y,d_totalforcex,d_totalforcey,th1);
	ReduceForces(d1_list1,d1_list2,d1_list3,d1_list4,d1_forces1x,d1_forces2x,d1_forces1y,d1_forces2y,d_totalforcex1,d_totalforcey1,th1);	
	 cout << "d1" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();	

	double *d_totalforcex2;
	double *d_totalforcey2;

	cudaMalloc((void**)&d_totalforcex2,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey2,totalN*sizeof(double));

	resetforce(d_totalforcex2);
	resetforce(d_totalforcey2);	 
	//ReduceForces(d2_list1,d2_list2,d2_list3,d2_list4,d2_forces1x,d2_forces2x,d2_forces1y,d2_forces2y,d_totalforcex,d_totalforcey,th2);
	ReduceForces(d2_list1,d2_list2,d2_list3,d2_list4,d2_forces1x,d2_forces2x,d2_forces1y,d2_forces2y,d_totalforcex2,d_totalforcey2,th2);
	cout << "d2" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex3;
	double *d_totalforcey3;

	cudaMalloc((void**)&d_totalforcex3,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey3,totalN*sizeof(double));

	resetforce(d_totalforcex3);
	resetforce(d_totalforcey3);	 	
	//ReduceForces(d3_list1,d3_list2,d3_list3,d3_list4,d3_forces1x,d3_forces2x,d3_forces1y,d3_forces2y,d_totalforcex,d_totalforcey,th3);	 	
	ReduceForces(d3_list1,d3_list2,d3_list3,d3_list4,d3_forces1x,d3_forces2x,d3_forces1y,d3_forces2y,d_totalforcex3,d_totalforcey3,th3);
	 cout << "d3" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex4;
	double *d_totalforcey4;

	cudaMalloc((void**)&d_totalforcex4,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey4,totalN*sizeof(double));

	resetforce(d_totalforcex4);
	resetforce(d_totalforcey4);	 		
	//ReduceForces(d4_list1,d4_list2,d4_list3,d4_list4,d4_forces1x,d4_forces2x,d4_forces1y,d4_forces2y,d_totalforcex,d_totalforcey,th4);
	ReduceForces(d4_list1,d4_list2,d4_list3,d4_list4,d4_forces1x,d4_forces2x,d4_forces1y,d4_forces2y,d_totalforcex4,d_totalforcey4,th4);
	 cout << "d4" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex5;
	double *d_totalforcey5;

	cudaMalloc((void**)&d_totalforcex5,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey5,totalN*sizeof(double));

	resetforce(d_totalforcex5);
	resetforce(d_totalforcey5);	 	 
	//ReduceForces(d5_list1,d5_list2,d5_list3,d5_list4,d5_forces1x,d5_forces2x,d5_forces1y,d5_forces2y,d_totalforcex,d_totalforcey,th5);
	ReduceForces(d5_list1,d5_list2,d5_list3,d5_list4,d5_forces1x,d5_forces2x,d5_forces1y,d5_forces2y,d_totalforcex5,d_totalforcey5,th5);	

	 cout << "d5" << endl;	
	// print_device_weave_float2(d5_list1,d5_list2,d5_forces1x,d5_forces1y,d_particles,th5,totalN);
	// pausel();
	double *d_totalforcex6;
	double *d_totalforcey6;

	cudaMalloc((void**)&d_totalforcex6,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey6,totalN*sizeof(double));

	resetforce(d_totalforcex6);
	resetforce(d_totalforcey6);	 	
	//ReduceForces(d6_list1,d6_list2,d6_list3,d6_list4,d6_forces1x,d6_forces2x,d6_forces1y,d6_forces2y,d_totalforcex,d_totalforcey,th6);	 	
	ReduceForces(d6_list1,d6_list2,d6_list3,d6_list4,d6_forces1x,d6_forces2x,d6_forces1y,d6_forces2y,d_totalforcex6,d_totalforcey6,th6);		
	 cout << "d6" << endl;
	// print_device_weave_float2(d6_list1,d6_list2,d6_forces1x,d6_forces1y,d_particles,th6,totalN);
	// cout << endl;
	// pausel();
	double *d_totalforcex7;
	double *d_totalforcey7;

	cudaMalloc((void**)&d_totalforcex7,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey7,totalN*sizeof(double));

	resetforce(d_totalforcex7);
	resetforce(d_totalforcey7);	 
	//ReduceForces3(d7_list1,d7_list2,d7_list3,d7_forces1x,d7_forces2x,d7_forces3x,d7_forces1y,d7_forces2y,d7_forces3y,d_totalforcex,d_totalforcey,th7);	 
	ReduceForces3(d7_list1,d7_list2,d7_list3,d7_forces1x,d7_forces2x,d7_forces3x,d7_forces1y,d7_forces2y,d7_forces3y,d_totalforcex7,d_totalforcey7,th7);	
	 cout << "d7" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex8;
	double *d_totalforcey8;

	cudaMalloc((void**)&d_totalforcex8,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey8,totalN*sizeof(double));
	

	
	resetforce(d_totalforcex8);
	resetforce(d_totalforcey8);	 
	//ReduceForces(d8_list1,d8_list2,d8_forces1x,d8_forces2x,d8_forces1y,d8_forces2y,d_totalforcex,d_totalforcey,th8);	 
	ReduceForces(d8_list1,d8_list2,d8_forces1x,d8_forces2x,d8_forces1y,d8_forces2y,d_totalforcex8,d_totalforcey8,th8);
	 cout << "d8" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double *d_totalforcex9;
	double *d_totalforcey9;

	cudaMalloc((void**)&d_totalforcex9,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey9,totalN*sizeof(double));

	resetforce(d_totalforcex9);
	resetforce(d_totalforcey9);	 
	//ReduceForces3(d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces2x,d9_forces3x,d9_forces1y,d9_forces2y,d9_forces3y,d_totalforcex,d_totalforcey,th9);		 
	ReduceForces3(d9_list1,d9_list2,d9_list3,d9_forces1x,d9_forces2x,d9_forces3x,d9_forces1y,d9_forces2y,d9_forces3y,d_totalforcex9,d_totalforcey9,th9);		
	
	double *d_totalforcex10;
	double *d_totalforcey10;

	cudaMalloc((void**)&d_totalforcex10,totalN*sizeof(double));
	cudaMalloc((void**)&d_totalforcey10,totalN*sizeof(double));

	resetforce(d_totalforcex10);
	resetforce(d_totalforcey10);	 
	 cout << "d9" << endl;	
	// print_device_array(d_totalforcex,totalN);
	// pausel();
	double ff = (v0_a+v0_b)/2.;
	//ReduceForcesAndNormalize(d10_list1,d10_forces1x,d10_forces1y,d_totalforcex,d_totalforcey, max_s, ff, th10);
	ReduceForcesAndNormalize(d10_list1,d10_forces1x,d10_forces1y,d_totalforcex10,d_totalforcey10, max_s, ff, th10);

	cout << "reduction" << endl;
	// cout << "d10" << endl;
	// print_device_array(d_totalforcex,totalN);
	// pausel();

	matrix<double> ftemp2(totalN,dimension),ftemp3(totalN,dimension);
	//matrix<double> angforces1(nc,dimension-1),angforces2(nc,dimension-1),angforces3(nc,dimension-1);
	// cout << "matrices initialized" << endl;

	matrix<double> F1((*obj).calculateforces(*froyo1,*faa)); //calculate the forces using the pairs as an input

	matrix<double> F2((*obj).calculateforces(*froyo2,*fbb)); //calculate the forces using the pairs as an input

	matrix<double> ftemp1((*obj).calculateforces(*froyo3,*fcc)); //calculate the forces using the pairs as an input
	
	matrix<double> F3((*obj).calculateforces(*froyo4,*fab)); //calculate the forces using the pairs as an input

	this->ForcesDueToPositionPL(*froyo5,ftemp2); //calculate the forces using the pairs as an input

	this->ForcesDueToPositionPL(*froyo6,ftemp3); //calculate the forces using the pairs as an input

	this->CalculateBindings(*froyo5,*froyo6);

	matrix<double> F4 = this->BindingForces();

	matrix<double> F5 = this->PositionForcesDueToAngles();


//	cout << "active forces" << endl;
	//cout << "pos forces" << endl;
	 matrix<double> F6((*obj).calculateforces(*bondpairs,*bindm));

	//cout << "bond forces" << endl;
	matrix<double> F7((*obj).calculateforces_threebody(*bendtriplets,*bendp));

	//cout << "after check matrix" << endl;

	matrix<double> F = ftemp1+ftemp2+ftemp3+F1+F2+F3+F4+F5+F6+F7;

	cout << l << endl;
	for(int j1 = 0 ; j1 < na+nb ; j1++ ) {
		if(bound[j1]>0) cout << j1 << ",";
	}
	cout << endl;
	print_device_array_indices(d_bound,na+nb);
	pausel();




	cout << F1 << endl;
	print_device_array_weave(d_totalforcex1,d_totalforcey1,totalN);
	cout << th1 << endl;
	cout << "faa" << endl;
	pausel();



	cout << F2 << endl;
	print_device_array_weave(d_totalforcex2,d_totalforcey2,totalN);
	cout << "fbb" << endl;
	pausel();


	cout << ftemp1 << endl;
	print_device_array_weave(d_totalforcex3,d_totalforcey3,totalN);
	cout << "fcc" << endl;
	pausel();		

	cout << F3 << endl;
	print_device_array_weave(d_totalforcex4,d_totalforcey4,totalN);
	cout << "fab" << endl;
	pausel();

	cout << ftemp2 << endl;
	print_device_array_weave(d_totalforcex5,d_totalforcey5,totalN);
	cout << "fac" << endl;
	pausel();		

	cout << ftemp3 << endl;
	print_device_array_weave(d_totalforcex6,d_totalforcey6,totalN);
	cout << "fbc" << endl;
	pausel();

	cout << F4 << endl;
	print_device_array_weave(d_totalforcex7,d_totalforcey7,totalN);
	cout << "bound to mt" << endl;
	pausel();	

	cout << F5 << endl;
	print_device_array_weave(d_totalforcex10,d_totalforcey10,totalN);
	cout << "position force" << endl;
	pausel();

	cout << F6 << endl;
	print_device_array_weave(d_totalforcex8,d_totalforcey8,totalN);
	cout << "bound within mt" << endl;
	pausel();		


	cout << F7 << endl;
	print_device_array_weave(d_totalforcex9,d_totalforcey9,totalN);
	cout << "bending force" << endl;
	pausel();		







	double *d_R1;
	double *d_R2;

	// cudaMalloc((void**)&d8_forces1x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces1y,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2x,th8*sizeof(double));
	// cudaMalloc((void**)&d8_forces2y,th8*sizeof(double));	
	cudaMalloc((void**)&d_R1,totalN*sizeof(double));
	cudaMalloc((void**)&d_R2,totalN*sizeof(double));

	setstaterandom(d_R1,1.732050808,totalN);
 	setstaterandom(d_R2,1.732050808,totalN);

 	print_device_array(d_R1,totalN);
 	print_device_array(d_R2,totalN);

 	pausel();

 	advmom2D(d_momenta, d_totalforcex, d_totalforcey, d_R1, d_R2, cons1,cons2,cons3,totalN);
 	advpos2D(d_particles, d_momenta, cons4, totalN);

 	applypbc2D(d_particles,d_momenta,l,is_periodic,totalN);

 	cout << "updated" << endl;

	cudaFree(d1_forces1x);
	cudaFree(d1_forces2x);
	cudaFree(d1_forces1y);
	cudaFree(d1_forces2y);
	cudaFree(d2_forces1x);
	cudaFree(d2_forces2x);
	cudaFree(d2_forces1y);
	cudaFree(d2_forces2y);
	cudaFree(d3_forces1x);
	cudaFree(d3_forces2x);
	cudaFree(d3_forces1y);
	cudaFree(d3_forces2y);
	cudaFree(d4_forces1x);
	cudaFree(d4_forces2x);
	cudaFree(d4_forces1y);
	cudaFree(d4_forces2y);
	cudaFree(d5_forces1x);
	cudaFree(d5_forces2x);
	cudaFree(d5_forces1y);
	cudaFree(d5_forces2y);
	cudaFree(d6_forces1x);
	cudaFree(d6_forces2x);
	cudaFree(d6_forces1y);
	cudaFree(d6_forces2y);
 	cudaFree(d7_list1);
 	cudaFree(d7_list2);
 	cudaFree(d7_list3);
	cudaFree(d7_forces1x);
	cudaFree(d7_forces2x);
	cudaFree(d7_forces3x);
	cudaFree(d7_forces1y);
	cudaFree(d7_forces2y);
	cudaFree(d7_forces3y);
	cudaFree(d8_forces1x);
	cudaFree(d8_forces2x);
	cudaFree(d8_forces1y);
	cudaFree(d8_forces2y);
	cudaFree(d9_forces1x);
	cudaFree(d9_forces2x);
	cudaFree(d9_forces3x);
	cudaFree(d9_forces1y);
	cudaFree(d9_forces2y);
	cudaFree(d9_forces3y);
	cudaFree(d10_list1);
	cudaFree(d10_forces1x);
	cudaFree(d10_forces1y);

	cout << "freed" << endl;


	}
}
*/