**Note**: LENT singularity has not been updated to work with OpenFOAM-v2112 yet.

## Using Singularity containers
1. First, build OpenFOAM version 1912 using the definition file 'openfoam-v1912.def' inside the 'singularity' folder.\
    `sudo singularity build openfoam.sif openfoam.def`
2. Then, build the container that consists of the environment and the build of the LENT project by using the definition file 'lent.def' inside the 'singularity' folder.\
    `sudo singularity build lent.sif lent.def`

### Private project LENT
If LENT is a private project, it cannot be accessed automatically. You need to modify the Definition file 'lent.def' in order to copy a pair of SSH-keys from your local
machine. It is advised that these keys be disposable, as they will continue to live inside the container. This modification will be removed once the project is public
and there will no need for authentication to clone.

### Use of the containers
Once the container is build you can execute the following commands on it:
1. clone, in order to get a copy of the project in your host system.\
    `./lent.sif clone`

2. build, in order to build the LENT project using the environment (dependencies) of the container.\
   `./lent.sif build lent`

Note: In the building command the last keyword 'lent' is refering to the name of the (host) project folder. It can be anything if you choose to rename that folder.
