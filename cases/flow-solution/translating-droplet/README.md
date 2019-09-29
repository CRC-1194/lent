# Translating droplet
Note: currently, only the three-dimensional variant is up-to-date and running

## Test background
Young-Laplace equilibrium for a spherical interface in the absence of external forces.
This case is very similar to the stationary droplet, except that the droplet is advected
with a constant, uniform background velocity. Thus, the interface advection plays a more important
role compared to the stationary droplet. However, the same exact solution is valid for this test case
in the moving reference frame of the droplet.  
This test case has first been proposed by Popinet:  
_Popinet, Stéphane. "An accurate adaptive solver for surface-tension-driven interfacial flows." Journal of Computational Physics 228.16 (2009): 5838-5866._  
Several setups from the literature are replicated (the corresponding publication
is referenced in the header of the parameter study). In addition, there are setups
using parameters from targeted "real-world" applications.

## What can be tested
* Well-balancedness of discretization/solution algorithm given an exact,
    constant curvature.
* Evaluation of capillary accuracy --> suitability of solution configuration
    for near-equilibrium applications.
* Test if the solution algorithm is able to reach/maintain a numerical equilibrium
    between surface tension and pressure.
* When compared to results from the stationary droplet: perturbation of equilibrium
    due to advection errors.

This test case is primarily a check for the interaction between  discretization, solution procedure
of the pressure-velocity coupling, numerical curvature model and interface advection.

## Error metrics
After subtraction of the background velocity field, the translating droplet test case is identical to
the stationary droplet. Thus, this test case also uses the _stationaryDropletFunctionObject_ from the _two-phase-validation_
library [TODO: insert hyperlink] to evaluate a set of error metrics. These are based on the
deviation of the numerically computed velocity field and pressure field from the
exact solution of the Young-Laplace equilibrium.

## Usage
Have a look at the README in the parent directory "flow-solution" to learn how to use
the hydrodynamic test cases.
