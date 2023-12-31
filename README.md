# LENT - a hybrid LEvel-set / froNT tracking (LENT) method in OpenFOAM 

The LENT OpenFOAM module implements the LENT hybrid unstructured Level Set / Front Tracking method for DNS simulations of two-phase flows using a collocated unstructured FVM equation discretization.

**IMPORTANT**: this is an actively developed research code. Report bugs [here](mailto:contact-project+leia-methods-lent-32701353-issue-@incoming.gitlab.com), contact us directly regarding  possible features or applications. 

## Developers

Tomislav Maric maric{a}mma[dot]tu-darmstadt[dot]de

Tobias Tolle tolle{a}mma[dot]tu-darmstadt[dot]de

Jun Liu liu{a}mma[dot]tu-darmstadt[dot]de

Mathematical Modeling and Analysis Institute,Mathematics Department, TU Darmstadt

## Publications

Liu, Jun, Tobias Tolle, Dieter Bothe, and Tomislav Maric. "A collocated unstructured finite volume Level Set/Front Tracking method for two-phase flows with large density-ratios." arXiv preprint arXiv:2109.01595 (2021). [arXiv:2109.01595](https://arxiv.org/abs/2109.01595)

Tolle, Tobias, Dieter Bothe, and Tomislav Marić. "SAAMPLE: A segregated accuracy-driven algorithm for multiphase pressure-linked equations." Computers & Fluids 200 (2020): 104450. [doi:10.1016/j.compfluid.2020.104450](https://doi.org/10.1016/j.compfluid.2020.104450) [	arXiv:2001.09775](https://arxiv.org/abs/2001.09775)

Marić, Tomislav, Holger Marschall, and Dieter Bothe. "lentFoam–A hybrid Level Set/Front Tracking method on unstructured meshes." Computers & Fluids 113 (2015): 20-31. [doi: 10.1016/j.compfluid.2014.12.019](https://doi.org/10.1016/j.compfluid.2014.12.019)

## Installation 

### Prerequisites

* OpenFOAM versions: [OpenFOAM v2112 (ESI)](https://develop.openfoam.com/Development/openfoam/-/blob/master/doc/Build.md)
* Compilers: gcc 10.3.0 and earlier
* Build system: CMake 3.14 or higher
* git for fetching dependencies
* gmsh or paraview for surface mesh generation
* admesh for normal consistency (only for some 2D test cases)

### Compilation and installation

LENT is using the [CMake](https://cmake.org) build system.
Make sure that you have sourced OpenFOAM's `bashrc` so that the OpenFOAM environment variables are set properly.
Then execute the following commands inside the `lent` directory:
    
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=./ -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=on ..
    make && make install

This will configure, compile and install the `lent` library
and its applications to `$FOAM_USER_LIBBIN` and `$FOAM_USER_APPBIN`, respectively.
Valid build types are `Release`, `Debug` and `RelWithDebInfo`. The flag
`-DCMAKE_EXPORT_COMPILE_COMMANDS=on` is optional.  
Alternatively, you can run the `install.sh` script given that `gcc-10` and `g++-10` are available on your
system.  
Lent comes with some scripts, e.g. to setup and run parameter studies. To be able to use them, you need to add the
`scripts` directory to your `PATH` variable by executing the following command in your shell:

    echo 'export PATH=$HOME/OpenFOAM/openfoam/lent/scripts:$PATH' >> $HOME/.bashrc
    . $HOME/.bashrc

This command assumes that Lent is located in `~/OpenFOAM/openfoam`. If this is not the case, change the path accordingly.

### Installing PyFoam 

PyFoam scripts are used directly by some LENT utilities, so PyFoam needs to be installed on the system. Information on PyFoam installation can be found on the [OpenFOAM Wiki](http://openfoamwiki.net/index.php/Contrib/PyFoam). 

## Usage

There are two important utilities that are distributed with the LENT repository which are used to prepare the simulation case:  

* `lentClearCasesRecursive` - `pyFoam` based script that cleans the LENT cases recursively
* `lentSetFields` - pre-processing application that sets the signed distance and cell search fields required by the LENT method
* `lentCreateFront` - create a surface STL file based on an analytical surface description.

for the lent solvers:

* `lentReconstruct` - solver application used to *reconstruct* the LENT front 
* `lentAdvect` - solver application used to *advect* the LENT front with a *prescribed velocity field*
* `lentFoam` - solver application used to execute two-phase Direct Numerical Simulations with the LENT method

Other testing and utility applications are available in the `applications` folder with the appropriate descriptions placed in the implementation file headers.  
In `cases/flow-solution/translating-droplet/three-dimensional` there are several run scripts providing examples
on the usage of scripts and solvers.

### Utilities information 

#### `lentClearCasesRecursive` 

* Expects to find `0.org` directory with the original fields in the simulation case directory. The LENT method modifies the initial `0` time-step directory, as it requires the signed and search distance fields to be computed before the start of the simulation. 
* Is used mostly for validation cases of the `lent-validation-cases` repository. 
* Even if `0.org` doesn't exist, it will clean up the `front/front-000*.vtk` files written from the last run. 

#### `lentSetFields` 

* Requires a file named `front.stl` in the folder `front` of the LENT simulation case. 
* It uses the surface mesh stored in the `front.stl` file to pre-process the signed and search distance fields. 

#### `lentCreateFront`

* Requires a dictionary named `frontSurface` in `system/lentSolution` specifying the surface.
* It uses the reconstruction procedure to generate a triangular surface mesh for the given surface
    and saves it as `constant/front.stl`.


### Running cases

This is the standard workflow: 

    lentClearCasesRecursive
    lentSetFields
    lentFoam 

`lentFoam` can be replaced by another solver application (`lentReconstruct` or `lentAdvect`).
