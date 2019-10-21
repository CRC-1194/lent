# Working with the rising bubble test case

Manually executing the case:

1. Initialize the study with "lentParameterInit -c templateCase -p parameterFile -s studyName".
2. In each parameter case, create the mesh with blockMesh.
3. In each parametr case, create the front with lentCreateFront.
4. In each parameter case, set the fields with lentSetFields.
5. In each parameter csae, use "topoSet". This creates the interface region based on the bubble position and radius to be refined with hex refinement. 
6. In each parameter case, refine the hex mesh by executing "refineHexMesh.sh". 

All this is done automatically by the runParameterStudy.sh in serial. On a cluster, a batch system can be used for the simultaneous execution in every parameter case.
