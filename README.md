# ASG_CFD

This CFD code has been developed by Alejandro Sanchez García as his final thesis for the Master's degree in Aerospace Engineering at Universitat Politècnica de Catalunya (UPC).
The final thesis is called *Study of the flow in a nozzle for rocket engines*, and it will be publicly available soon on [UPCommons](https://upcommons.upc.edu/).

This is a MATLAB code designed to be able to implement a computationally affordable, one-dimensional, compressible model for the flow of a rocket engine nozzle. To that end, three
numerical schemes have been implemented: the MacCormack technique, the Godunov method, and the MUSCL-Hancock scheme.

The purpose of developing such a code is to be able to integrate it into algorithms designed to optimize nozzle and supersonic air inlets, among other geometries, described by a certain
area distribution function $A = A(x)$. This may offer a preliminary design tool to derive unconventional shapes intended to enhance the performance of rockets and other flying vehicles.

## Main functions

The main functions that constitute the numerical schemes are the following.

### Solvers

- *SOLVER_MacCormack_time*: implements the time-marching MacCormack technique.
- *SOLVER_MacCormack_space*: implements the space-marching MacCormack technique.
- *SOLVER_Godunov*: implements the second version of Godunov's first-order upwind method.
- *SOLVER_MUSCL_Hancock*: implements the second-order MUSCL-Hancock scheme.

Except for the *SOLVER_MacCormack_space*, the user can impose a maximum number of time steps for all schemes.

### Mesh generation

- *generate_mesh*: discretises the computational domain according to a certain number of nodes $N$, to be given by the user.

### Steady state convergence

- *check_convergence*: computes the residual of the whole density field between current and last time steps. Once a certain tolerance $TOL$ is met (can be modified by the user), the steady state is reached.

### Slope vectors

- *slope_vector*: computes the slope vectors required by the MUSCL-Hancock scheme. The user can select a TVD method to be used, among the following:
  - *none*
  - *limited_slopes*
  - *slope_limiter_SUPERBEE*
  - *slope_limiter_MINBEE*
 
## Auxiliary functions

This involves all the remaining scripts, which include functions defining the state equation being used (ideal gas law), functions to get the gas properties ($R_g$, $\gamma$ and $c_v$), functions that implement
the exact and the HLLC approximate solvers to the Riemann problem, among others.

Further information can be found in the *Annexes* to the thesis, as well as in the functions own documentation by typing *help function_name* in MATLAB.

## Validation scripts

Folder **Test_cases** includes the main scripts describing the different problems that have been used to test and validate the code: the Sod shock tube problem, a Prandtl-Meyer expansion, and different nozzle flow cases (a subsonic-supersonic expansion, a case with
an inner shock wave, and a purely subsonic expansion). These nozzle cases are reported by Anderson [1], whose numerical results have been tabulated in three different text files.

Note that script *MAIN_Nozzle_HGS* may require [HGS](https://github.com/ManelSoria/HGS), an open-source thermochemistry calculator for mixtures of ideal gases at high temperature, which has been integrated into the CFD code to implement two gas models in addition
to the **perfect gas** model: **frozen flow** and **shifting flow**.

Again, further information can be found in the *Annexes* to the thesis.

## References

[1]  ANDERSON, J. D. *Computational Fluid Dynamics: The Basics with Applications*. New York: McGraw-Hill, 1995. ISBN 0070016852.
