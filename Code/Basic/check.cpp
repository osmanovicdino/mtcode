

	/*
	savef1forces.push_back(F1);
	savef2forces.push_back(F2);
	savef3forces.push_back(F3);
	savef4forces.push_back(F4);
	savef5forces.push_back(F5);
	savef6forces.push_back(F6);
	savef7forces.push_back(F7);
	saveftemp1forces.push_back(ftemp1);
	saveftemp2forces.push_back(ftemp2);
	saveftemp3forces.push_back(ftemp3);
	matrix<double> temppos = obj->getdat();
	savepositions.push_back(temppos);
	savebound.push_back(bound);
	saveboundalongs.push_back(bound_along);
	int mxsiz =  25;

	if(savef1forces.size() > mxsiz ) {
		savepositions.erase(savepositions.begin());
		savebound.erase(savebound.begin());
		saveboundalongs.erase(saveboundalongs.begin());
		savef1forces.erase(savef1forces.begin());
		savef2forces.erase(savef2forces.begin());
		savef3forces.erase(savef3forces.begin());
		savef4forces.erase(savef4forces.begin());
		savef5forces.erase(savef5forces.begin());
		savef6forces.erase(savef6forces.begin());
		savef7forces.erase(savef7forces.begin());
		saveftemp1forces.erase(saveftemp1forces.begin());
		saveftemp2forces.erase(saveftemp2forces.begin());
		saveftemp3forces.erase(saveftemp3forces.begin());
	}
	*/




	/*
	double maxsize = 4000;
	if(chckmatrixsize(F1,maxsize)) {
		cout << "F1" << endl;
		chckmatrixprint(F1,maxsize);
		vector1<int> indices =  chckmatrixindices(F1,maxsize);
		//cout << indices << endl;
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}		
		error("error in F1");
	}
	if(chckmatrixsize(F2,maxsize)) {
		cout << "F2" << endl;
		chckmatrixprint(F2,maxsize);
		vector1<int> indices =  chckmatrixindices(F2,maxsize);
		//cout << indices << endl;
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}		
		error("error in F2");
	}
	if(chckmatrixsize(F3,maxsize)) {
		cout << "F3" << endl;
		chckmatrixprint(F3,maxsize);
		vector1<int> indices =  chckmatrixindices(F3,maxsize);
		//cout << indices << endl;
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}		
		error("error in F3");
	}
	if(chckmatrixsize(F4,maxsize)) {
		cout << "F4" << endl;
		chckmatrixprint(F4,maxsize);
		//chckmatrixprint(savebindingforces.at(savebindingforces.size()-2),100);
		vector1<int> indices =  chckmatrixindices(F4,maxsize);
		//cout << indices << endl;
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}

		cout << "old binding forces" << endl;
		this->BindingForcesTest(savebound[savebound.size()-2],saveboundalongs[saveboundalongs.size()-2],savepositions[savepositions.size()-2]);


		cout << "new binding forces" << endl;
		this->BindingForcesTest(savebound[savebound.size()-1],saveboundalongs[saveboundalongs.size()-1],savepositions[savepositions.size()-1]);

		error("error in F4");
	}
	if(chckmatrixsize(F5,maxsize)) {
		cout << "F5" << endl;
		chckmatrixprint(F5,maxsize);
		vector1<int> indices =  chckmatrixindices(F5,maxsize);
		//cout << indices << endl;
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}		
		error("error in F5");
	}
	if(chckmatrixsize(F6,maxsize)) {
		cout << "F6" << endl;
		chckmatrixprint(F6,maxsize);
		vector1<int> indices =  chckmatrixindices(F6,maxsize);
		//cout << indices << endl;
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}		
		error("error in F6");
	}
	if(chckmatrixsize(F7,maxsize)) {
		cout << "F7" << endl;
		chckmatrixprint(F7,maxsize);
		vector1<int> indices =  chckmatrixindices(F7,maxsize);
		//cout << indices << endl;
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}		
		error("error in F7");
	}
	if(chckmatrixsize(ftemp1,maxsize)) {
		cout << "ftemp1" << endl;
		chckmatrixprint(ftemp1,1000);
		//chckmatrixprint(saveflforces1.at(saveflforces1.size()-2),1000);
		vector1<int> indices =  chckmatrixindices(ftemp1,1000);
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}
		error("error in ftemp2");

	}	

	if(chckmatrixsize(ftemp2,maxsize)) {
		cout << "ftemp2" << endl;
		chckmatrixprint(ftemp2,1000);
		//chckmatrixprint(saveflforces1.at(saveflforces1.size()-2),1000);
		vector1<int> indices =  chckmatrixindices(ftemp2,1000);
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}
		error("error in ftemp2");

	}

	if(chckmatrixsize(ftemp3,maxsize)) {
		cout << "ftemp3" << endl;
		chckmatrixprint(ftemp3,1000);
		//chckmatrixprint(saveflforces2.at(saveflforces2.size()-2),1000);
		vector1<int> indices =  chckmatrixindices(ftemp3,1000);
		for(int j = 0  ; j < savef1forces.size() ; j++) {
			cout << "positions: " << j << " " <<  savepositions[j](indices[0],0) << " " << savepositions[j](indices[0],1) << " " << savepositions[j](indices[1],0) << " " << savepositions[j](indices[1],1) << endl;
			cout << "distance: "  << (obj->getgeo()).distance(savepositions[j],indices[0],indices[1]) << endl;
			cout << "forces1: " << j << " " << savef1forces[j](indices[0],0) << " " << savef1forces[j](indices[0],1) << " " << savef1forces[j](indices[1],0) << " " << savef1forces[j](indices[1],1) << endl;
			cout << "forces2: " << j << " " << savef2forces[j](indices[0],0) << " " << savef2forces[j](indices[0],1) << " " << savef2forces[j](indices[1],0) << " " << savef2forces[j](indices[1],1) << endl;
			cout << "forces3: " << j << " " << savef3forces[j](indices[0],0) << " " << savef3forces[j](indices[0],1) << " " << savef3forces[j](indices[1],0) << " " << savef3forces[j](indices[1],1) << endl;
			cout << "forces4: " << j << " " << savef4forces[j](indices[0],0) << " " << savef4forces[j](indices[0],1) << " " << savef4forces[j](indices[1],0) << " " << savef4forces[j](indices[1],1) << endl;
			cout << "forces5: " << j << " " << savef5forces[j](indices[0],0) << " " << savef5forces[j](indices[0],1) << " " << savef5forces[j](indices[1],0) << " " << savef5forces[j](indices[1],1) << endl;
			cout << "forces6: " << j << " " << savef6forces[j](indices[0],0) << " " << savef6forces[j](indices[0],1) << " " << savef6forces[j](indices[1],0) << " " << savef6forces[j](indices[1],1) << endl;
			cout << "forces7: " << j << " " << savef7forces[j](indices[0],0) << " " << savef7forces[j](indices[0],1) << " " << savef7forces[j](indices[1],0) << " " << savef7forces[j](indices[1],1) << endl;
			cout << "forcestemp1: " << j << " " << saveftemp1forces[j](indices[0],0) << " " << saveftemp1forces[j](indices[0],1) << " " << saveftemp1forces[j](indices[1],0) << " " << saveftemp1forces[j](indices[1],1) << endl;
			cout << "forcestemp2: " << j << " " << saveftemp2forces[j](indices[0],0) << " " << saveftemp2forces[j](indices[0],1) << " " << saveftemp2forces[j](indices[1],0) << " " << saveftemp2forces[j](indices[1],1) << endl;
			cout << "forcestemp3: " << j << " " << saveftemp3forces[j](indices[0],0) << " " << saveftemp3forces[j](indices[0],1) << " " << saveftemp3forces[j](indices[1],0) << " " << saveftemp3forces[j](indices[1],1) << endl;
			cout << "binding: " << savebound[j][indices[0]] << " " << saveboundalongs[j][indices[0]] << endl;
			cout << endl;
		}
		error("error in ftemp3");
	}
	*/

/*
	if(chckmatrix(F1)) {
		error("F1");
	}
	if(chckmatrix(F2)) {
		error("F2");
	}
	if(chckmatrix(F3)) {
		error("F3");
	}
	if(chckmatrix(F4)) {
		error("F4");
	}

	if(chckmatrix(F6)) {
		error("F6");
	}
	if(chckmatrix(F7)) {
		error("F7");
	}
	if(chckmatrix(F5)) {
		error("F5");
	}
	if(chckmatrix(ftemp1)) {
		error("ftemp1");
	}
	if(chckmatrix(ftemp2)) {
		error("ftemp2");
	}
	if(chckmatrix(ftemp3)) {
		error("ftemp3");
	}
*/