# LENT - a hybrid LEvel-set / froNT tracking (LENT) method in OpenFOAM 

This is the library and application toolbox that implements the LENT method for DNS simulations of two-phase flows. 

Developers 

Tomislav Maric maric{a}mma[dot]tu-darmstadt[dot]de

Tobias Tolle tolle{a}mma[dot]tu-darmstadt[dot]de

Mathematics Department, Mathematical Modeling and Analysis Institute, TU Darmstadt

## Installation 

### Prerequisites

* OpenFOAM versions: [OpenFOAM-plus v1912](https://openfoam.com/releases/openfoam-v1912/)
* Compilers: Gcc 8 and above (tested with 9.2.1)
* Build system: CMake 3.14 or higher
* Git for fetching dependencies
* gmsh or paraview for surface mesh generation
* admesh for normal consistency 

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

for the lent solvers:

* `lentReconstruct` - solver application used to *reconstruct* the LENT front 
* `lentAdvect` - solver application used to *advect* the LENT front with a *prescribed velocity field*
* `lentFoam` - solver application used to execute two-phase Direct Numerical Simulations with the LENT method

Other testing and utility applications are available in the `applications` folder with the appropriate descriptions placed in the implementation file headers. 

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

