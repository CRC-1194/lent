# LENT - a hybrid LEvel-set / froNT tracking (LENT) method in OpenFOAM 

This is the library and application toolbox that implements the LENT method for DNS simulations of two-phase flows. 

Developers 

Tomislav Maric maric{a}mma[dot]tu-darmstadt[dot]de

Tobias Tolle tolle{a}mma[dot]tu-darmstadt[dot]de

Mathematics Department, Mathematical Modeling and Analysis Institute, TU Darmstadt

## Installation 

### Prerequisites

* OpenFOAM versions: 2.2.x  
* Compilers: gcc-4.8.2 and gcc-4.9.0
* `subversion` and `cmake` for compiling Google Test 

### Compiling Google Test

* Make sure that you have configured the OpenFOAM environment:

    ls $WM_PROJECT_DIR/etc/bashrc

should output the location of the `bashrc` OpenFOAM configuration file for your chosen version.  

* Source the *bashrc* configuration script from within the `lent` directory:

    . etc/bashrc

*Note* - CXXFLAGS variable is taken over from OpenFOAM by this step, so if skipped, GTest will be compiled with different compiler flags than LENT which may lead to problems.

* Prepare the Google Test third party directory 

    mkdir $LENT_PROJECT/third-party

* Get the Google Test sources 

    cd $LENT_PROJECT/third-party

    svn checkout http://googletest.googlecode.com/svn/trunk/ gtest

* Build Google Test with cmake 

    cd gtest

    mkdir build

    cd build 

    cmake -DBUILD_SHARED_LIBS=ON -Dgtest_build_samples=ON -G"Unix Makefiles" ..

    make


*Note* - if everything above was successful, LENT library, solvers and testing applications can be compiled. If not, you should still be able to compile the library and the solvers, and the testing applications will report errors during compilation.

### Compiling the LENT library and solvers 

*Note* : the following installation instructions do not compile the testing applications: if Google Test has not been built alongside LENT as described above, testing applications will fail to compile. 

* Source the *bashrc* configuration script from within the top directory:

    . etc/bashrc

* Execute Allwmake from within the top directory:

    . ./Allwmake

For permanent configuration of the environmental variables, source the *bashrc* script from your local *~/.bashrc*: 

    source path/to/lent/etc/bashrc

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


### Running cases

This is the standard workflow: 

    lentClearCasesRecursive
    lentSetFields
    lentFoam 

`lentFoam` can be replaced by another solver application (`lentReconstruct` or `lentAdvect`).

