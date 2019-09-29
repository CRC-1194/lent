# Stationary droplet
Note: currently, only the three-dimensional variant is up-to-date and running

## Test background
Young-Laplace equilibrium for a spherical interface in the absence of external forces.
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

This test case primarily a check for the discretization, the solution procedure
of the pressure-velocity coupling and the numerical curvature model.

## Error metrics
This test case uses the _stationaryDropletFunctionObject_ from the _two-phase-validation_
library [TODO: insert hyperlink] to evaluate a set of error metrics. These are based on the
deviation of the numerically computed velocity field and pressure field from the
exact solution of the Young-Laplace equilibrium.

## Usage
Have a look at the README in the parent directory "flow-solution" to learn how to use
the hydrodynamic test cases.
