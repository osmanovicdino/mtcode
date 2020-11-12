#ifndef LANGEVINRFORCE_CPP
#define LANGEVINRFORCE_CPP

void LangevinNVTR::calculate_forces_and_torques3D_onlyone(matrix<int> &pairs, ComboPatch &iny, BinaryBindStore &bo, AbstractBindingModel &bm, matrix<double> &forces, matrix<double> &torques)
{

    //for all the pairs, for all bindings

    //for a given sphere geometry

    //int np1 = sqrt(iny.getsize());
    int total_number_of_patches = iny.get_total_patches(this->getN());


    vector1<int> tempbound(total_number_of_patches); //no binding to begin wtih

    int depth_of_matrix = 4; //Choose this value to be deep enough such that all values can be stored

    matrix<int> boindices(total_number_of_patches, depth_of_matrix);

    // if (SQR(np1) != iny.getsize())
    // {
    //     error("number of patches and size of potential bundle incorrect");
    // }

    //vector1<int> nump(this->getN(),0); //this is the number of bound particles per particle

#pragma omp parallel for
    for (int i = 0; i < pairs.getNsafe(); ++i)
    {
        int p1 = pairs(i, 0);
        int p2 = pairs(i, 1);
        //int i1 = pairs(i,2);
        double dis;
        //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
        vector1<double> un(dimension);
        geo->distance_vector(*dat, p1, p2, un, dis);

        //un = i-j
        dis = sqrt(dis);

        un /= dis;

        double qtemp0 = orient->gpcons(p1, 0);
        double qtemp1 = orient->gpcons(p1, 1);
        double qtemp2 = orient->gpcons(p1, 2);
        double qtemp3 = orient->gpcons(p1, 3);
        double qtemp4 = orient->gpcons(p1, 4);
        double qtemp5 = orient->gpcons(p1, 5);
        double qtemp6 = orient->gpcons(p1, 6);
        double qtemp7 = orient->gpcons(p1, 7);
        double qtemp8 = orient->gpcons(p1, 8);

        double gtemp0 = orient->gpcons(p2, 0);
        double gtemp1 = orient->gpcons(p2, 1);
        double gtemp2 = orient->gpcons(p2, 2);
        double gtemp3 = orient->gpcons(p2, 3);
        double gtemp4 = orient->gpcons(p2, 4);
        double gtemp5 = orient->gpcons(p2, 5);
        double gtemp6 = orient->gpcons(p2, 6);
        double gtemp7 = orient->gpcons(p2, 7);
        double gtemp8 = orient->gpcons(p2, 8);

        // for (int j = 0; j < iny.num_patches(p1) ; j++)
        // {
        //     for (int k = 0; k < iny.num_patches(p2); k++)
        //     {

                //int potn = np1 * j + k;
            iny.UpdateIterator(p1,p2);
            for (int tp = 1; tp < (*(iny.p))[0] + 1; tp++)
            {
                int potn = (*(iny.p))[tp];
                // vector1<double> params = (iny.potential_bundle)[potn]->getparameters();

                int nxb1;// = params[0]; //iny[potn]->nxb1;
                int nyb1;// = params[1]; //iny[potn]->nyb1;
                int nzb1;// = params[2]; //iny[potn]->nzb1;

                int nxb2;// = params[3]; //iny[potn]->nxb2;
                int nyb2;// = params[4]; //iny[potn]->nyb2;
                int nzb2;// = params[5]; //iny[potn]->nzb2;

                double disp;// = params[6];

                double thetam;// = params[8];

                iny.get_params(p1,p2,potn,nxb1,nyb1,nzb1,nxb2,nyb2,nzb2,disp,thetam); //for this potential, get all the parameters

                double nx1 = nxb1 * qtemp0 + nyb1 * qtemp3 + nzb1 * qtemp6;
                double ny1 = nxb1 * qtemp1 + nyb1 * qtemp4 + nzb1 * qtemp7;
                double nz1 = nxb1 * qtemp2 + nyb1 * qtemp5 + nzb1 * qtemp8;

                double nx2 = nxb2 * gtemp0 + nyb2 * gtemp3 + nzb2 * gtemp6;
                double ny2 = nxb2 * gtemp1 + nyb2 * gtemp4 + nzb2 * gtemp7;
                double nz2 = nxb2 * gtemp2 + nyb2 * gtemp5 + nzb2 * gtemp8;

                double argthetai = -(nx1 * un.gpcons(0) + ny1 * un.gpcons(1) + nz1 * un.gpcons(2));
                double argthetaj = (nx2 * un.gpcons(0) + ny2 * un.gpcons(1) + nz2 * un.gpcons(2));

                if (argthetai > cos(thetam) && argthetaj > cos(thetam) && dis < 1.2 * disp)
                {
                    int wp1,wp2;
                    iny.which_patch(p1,p2,potn,wp1,wp2);
                    
                    boindices(wp1, tempbound[wp1]) = wp2;
                    boindices(wp2, tempbound[wp2]) = wp1;

                    tempbound[wp1] += 1;
                    tempbound[wp2] += 1;
                }
            }
        //     }
        // }
    }


    //Calculate all the bindings for all the patches

    vector1<int> indexes(total_number_of_patches);

    vector1<int> nbins = ConnectedComponents(boindices, tempbound, indexes);

    //Now we have the clusters. For each of these clusters, there is so some transition rate from one to another
    for (int i = 0; i < nbins.getsize() - 1; i++)
    {

        int size_of_cluster = nbins[i + 1] - nbins[i];

        if (size_of_cluster == 1)
        {
            //do nothing
            int i1 = indexes[nbins[i]];
            //isbound[i1]=false;

            bo.isbound[i1] = false;

            //not bound to anything.
        }
        else if (size_of_cluster == 2)
        {
            //all fine, bindings
            int ti1 = indexes[nbins[i]];
            int ti2 = indexes[nbins[i] + 1];

            int i1;
            int i2;
            sort_doublet(ti1,ti2,i1,i2);

            bool alreadybound_to_eachother = bo.boundto[i1] == i2 && bo.isbound[i1] && bo.isbound[i2];

            bool aft;

            bm.doublet(alreadybound_to_eachother, i1, i2, aft);

            if (aft)
            {
                bo.boundto[i1] = i2;
                bo.boundto[i2] = i1;
                bo.isbound[i1] = true;
                bo.isbound[i2] = true;
            }
            else
            {
                bo.isbound[i1] = false;
                bo.isbound[i2] = false;
            }
            //bool already_bound = prebound(i1,i2);
        }
        else if (size_of_cluster == 3)
        {
            //is the cluster fully connected or not?

            int ti1 = indexes[nbins[i]];
            int ti2 = indexes[nbins[i] + 1];
            int ti3 = indexes[nbins[i] + 2];

            //SORT THE INDICES (IMPORTANT)

            int i1;
            int i2;
            int i3;

            sort_triplet(ti1, ti2, ti3, i1, i2, i3);

            //bool b12,b23,b13;

            //DETERMINE WHETHER THEY ARE BOUND
            bool b12 = bo.boundto[i1] == i2 && bo.isbound[i1] && bo.isbound[i2];
            bool b23 = bo.boundto[i2] == i3 && bo.isbound[i2] && bo.isbound[i3];
            bool b13 = bo.boundto[i1] == i3 && bo.isbound[i1] && bo.isbound[i3];

            //DETERMINE THE CONNECTIVENESS OF THE GRAPH
            //remember, that in order to count as a triplet
            bool c12 = false;
            bool c23 = false;
            bool c13 = false;

            int nb1 = tempbound[i1];
            if (nb1 == 1)
            {
                int tempi = boindices(i1, 0);
                if (tempi == i2)
                {
                    c12 = true;
                    c13 = false;
                    c23 = true; //in order to be a triplet
                }
                else if (tempi == i3)
                {
                    c13 = true;
                    c12 = false;
                    c23 = true;
                }
                else
                    error("something weird");

                //check the other
            }
            else if (nb1 == 2)
            {
                c12 = true;
                c13 = true;

                int nb2 = tempbound[i2];
                if (nb1 == 1)
                {
                    c23 = false;
                }
                else
                {
                    c23 = true;
                }
            }
            else
            {
                error("error in clustering algorithm");
            }

            bool a12;
            bool a23;
            bool a13;
            bm.triplet(b12, b23, b13, c12, c23, c13, i1, i2, i3, a12, a23, a13);

            if (a12)
            {
                bo.boundto[i1] = i2;
                bo.boundto[i2] = i1;
                bo.isbound[i1] = true;
                bo.isbound[i2] = true;
                bo.isbound[i3] = false;
            }
            else if (a23)
            {
                bo.boundto[i2] = i3;
                bo.boundto[i3] = i2;
                bo.isbound[i1] = false;
                bo.isbound[i2] = true;
                bo.isbound[i3] = true;
            }
            else if (a13)
            {
                bo.boundto[i1] = i3;
                bo.boundto[i3] = i1;
                bo.isbound[i1] = true;
                bo.isbound[i2] = false;
                bo.isbound[i3] = true;
            }
            else
            {
                bo.isbound[i1] = false;
                bo.isbound[i2] = false;
                bo.isbound[i3] = false;
            }
        }
        else
        {
            error("the mythical 4 cluster, get to work, slob!");
        }
    }

    //Having gone through all the patches to determine the bindings, we can now calculate the forces

    vector1<bool> visited(total_number_of_patches); //keep track of patches that we already calculated

    for (int i = 0; i < total_number_of_patches ; ++i)
    {
        if (bo.isbound[i] == true && visited[i] == false) //only for bound patches we haven't visisted do we calculate forces
        {
            visited[i] = true; // We have now visisted this patch
            visited[bo.boundto[i]] = true; //we have also visisted the patch that it is bound to
            
            // int p1 = floor(i / np1);             //particle number 1
            // int p2 = floor(bo.boundto[i] / np1); //particle number 2

            //get the particle numbers from the patch numbers:

            int p1;
            int p2;

            iny.which_particle(i,bo.boundto[i],p1,p2);

            double dis;
            //vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
            vector1<double> un(dimension);
            geo->distance_vector(*dat, p1, p2, un, dis);

            //un = i-j

            int potn = iny.which_potential(p1,p2,i,bo.boundto[i]);
           // int potn = (i % np1) * np1 + (bo.boundto[i] % np1);
            dis = sqrt(dis);

            un /= dis;

            double fx;
            double fy;
            double fz;

            double tix;
            double tiy;
            double tiz;

            double tjx;
            double tjy;
            double tjz;

            (iny.potential_bundle)[potn]->force_and_torque(un, dis, *orient, p1, p2, fx, fy, fz, tix, tiy, tiz, tjx, tjy, tjz);

            forces(p1, 0) += fx;
            forces(p1, 1) += fy;
            forces(p1, 2) += fz;

            forces(p2, 0) += -fx;
            forces(p2, 1) += -fy;
            forces(p2, 2) += -fz;

            torques(p1, 0) += tix;
            torques(p1, 1) += tiy;
            torques(p1, 2) += tiz;

            torques(p2, 0) += tjx; // - dis * (fz * un[1] - fy * un[2]);
            torques(p2, 1) += tjy; // - dis * (fz * un[0] + fx * un[2]);
            torques(p2, 2) += tjz; // - dis * (fy * un[0] - fx * un[1]);
        }
        else
        {
            //do nothing
        }
    }

    //now we have only the real forces, we no longer need to calculate the forces for the non-bound particles:
}

#endif /* LANGEVINRFORCE_CPP */
