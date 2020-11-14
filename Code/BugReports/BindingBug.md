Observing a bug with binding when the two patches have the same body orientation.

I.e. when nx=1,ny=0,nz=0 for both particles, the torques appear to destroy the possibility of binding. (at finite temperature)

When nx1=1.0, ny1=0.0, nz1 =0.0 and nx2 =-1.0 and ny2= 0.0 and nz2 =0.0 this bug does not exist (at finite temp)

for starting conditions where p1 = (0,0,0) and p2 = (x,0,0) for x less than d (where d is the binding)

they are oriented to face each other (orientations identity matrix for both)

this case works ok.


for the nx1 = 1.0, ny1 =0.0, nz1 = 0.0 and nx2 = 1.0, ny=0.0, nz =0.0

with orientation 1 an identity matrix and orientation2 is (-1,0,0,0,1,0,0,0,-1) (pointing back) we get this effect.

Attempt 1:

for the different bond angles, try interchanging positions and orienting them negatively.

THIS FAILS

Now try orienting them both positively, but swapping nx1 and nx2 in the potential.

THIS WORKS

Now try to cases where it fails, but take un->-un

THIS DOES NOT WORK, PARTICLES DO NOT SEE EACH OTHER, DUE TO ANGLE CONDITION

Try with other but with un after angle condition

THIS LEADS TO REPULSION OF FX

Try for only the torques:

THIS LEADS TO INCORRECT TORQUES, BUT NOW WITH ROTATION IN Z AND NOT IN Y

Check again the torque calculations for your potentials. I.e. track the torque and the force at the same time.

Systematic error in the lab frame torques (it appears). Going to try reducing the viscosity.

This makes effect stronger

Going to try at kt=0.0. For a perfectly aligned sphere it works. Try off alignment,

for 
<pre><code>
    double nxtemp = 0.95;
    double nytemp = 0.31225;
    KernFrenkelOnePatch2 testpot(nxtemp, nytemp, 0., -nxtemp, -nytemp, 0., 100., 2., pi / 3., 0.75);
</code></pre>
   
the entire system rotates
<pre><code>
    double nxtemp = 0.95;
    double nytemp = 0.31225;
    KernFrenkelOnePatch2 testpot(nxtemp, nytemp, 0., -nxtemp, nytemp, 0., 100., 2., pi / 3., 0.75);
</code></pre>
The patches are repelled but do not rotate

all the cases where the nz are off behave correctly, problem with ny

when thw two ny are in the same direction, the torques are wrong. When they are in a different direction the torques appear correct but the whole system rotates.

The issue appears to be in the order of the mutliplication of the generator of rotations and the orientations.

Changing from $Q(t) = Q(0) R^T$ to $Q(t) = R^T Q(0)$ appears to work. 

To determine this, I went back to basics with regards to the rotate() function.

$$ \frac{d \mathbf{u}(t)}{d t} = \mathbf{\omega}(t) \times \mathbf{u}(t) $$

Which can be rewritten in matrix form:

$$ \frac{d u(t)}{d t} = \underline{\Omega(t)}.u(t) $$
where we have now included $\Omega$ as a matrix. The soluton to this equation is given by:

$$\mathbf{u}(t)=\exp\left(\underline{\int\Omega} (t)\mathrm{d}t\right).\mathbf{u}(0)$$
this gives a matrix rule for updating the orientation in time (we can transform $\mathbf{u}$ into a matrix which represents coordinate transformations from lab fixeed frame to body fixed frame)

If the code fails or is bugged again, replace their complicated rotate formalism with this method.



