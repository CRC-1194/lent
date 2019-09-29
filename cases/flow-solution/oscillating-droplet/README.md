# Oscillating droplet
Note: currently, only the three-dimensional variant is up-to-date and running

## Test background
[TODO: do a more thorough literature survey and have a look at the book by Lamb].
For small perturbations of a spherical interface, there exist analytical solutions
by Lamb [TODO: add citation] for the oscillation frequency and the decay of the
amplitude (in case the viscosity is non-zero). So, numerical setups fulfilling
the small perturbations requirement can be compared quantitatively with these analytical solutions.  
However, qualitative conclusions can also be drawn by monitoring the evolution of surface energy
and kinetic energy. For inviscid fluids, these metrics show how well (or bad) a two-phase
method fulfills the conservation of energy.  
Due to the interface movement, this test allows to evaluate the mass conservation properties
of a method.

## What can be tested
* For inviscid fluids, test energy conservation of a method.
* Quantitative comparison with Lamb's analytical solutions for some setups.
* Conservation of mass or volume.
* Conservativeness of the surface tension model: does the droplet's centre of
    gravity move in the course of simulation?

Apart from breakup/coalescence this test case involves all parts of a two-phase method.


## Error metrics
This test case uses the _oscillatingDropletFunctionObject_ from the _two-phase-validation_
library [TODO: insert hyperlink] to evaluate a set of metrics. These are not error metrics
themselves, but energy metrics, the droplet's centre of gravity and its semi-axes.

## Usage
Have a look at the README in the parent directory "flow-solution" to learn how to use
the hydrodynamic test cases.
