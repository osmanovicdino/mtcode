In combopatch, we define methods for each particle being "patchy". They each have patches which can interact with one another on their surfaces.

As an example, consider a patchy particle with 4 patches, the 4 patches would lead to the 4*4=16 different potentials interacting between these two particles. We wish to store these potentials, and thus keep track of the dynamics on the level of the patches, the particle identites and the potentials.

In combopatch.h the abstract base class here is defined. The basic goal of this class is to:

-provide a container to store all the potentials in the problem. This will be a list of potentials.

-when two particles meet, select from the list of all the potentials the indices corresponding to the two particles. For example if my system has a particle with 4 patches and a particle with 2 patches we need a scheme whereby we select the correct potentials to loop over depnding on the identities of the particles. We store this as a list of indices that return the correct potentials from the container mentioned above. 

-For two particles of index i and j interacting with a potential index, potn, return the patches on the particles this corresponds to.

- For two patches, return the particles

- For two patches, return the potential index

In general, these last three problems will depend on how the potentials are stored. We use the following scheme for storing 4 patch particles, for instance:
(patch no on particle1, patch number on particle 2),(potential number)

(0,0) 0
(0,1) 1
(0,2) 2
(0,3) 3
(1,0) 4
(1,1) 5
(1,2) 6
(1,3) 7
(2,0) 8
(2,1) 9
(2,2) 10
(2,3) 11
(3,0) 12
(3,1) 13
(3,2) 14
(3,3) 15

so for example, if we are given patch nos (i,j) we can get the potential via (i*4)+j

what if we are instead given the potential number? The mapping is still unique:

i = floor(potn/4)
j = mod(potn,4)

this can be worked out for any combination of particles with different numbers of patches. This means that there is also a symmetry in patches. For example, let us say we have particle i and particle j with potential 12, this must be the same as particle j and particle i with potential 3