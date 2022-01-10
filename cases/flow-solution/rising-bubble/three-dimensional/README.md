# Working with the rising bubble test case

## Manually executing the case:

Enable the use of case scripts

```
three-dimensional > . ../../../scripts/bashrc
```

Create the parameter study for the Duineveld parameters using `blockMesh`:

```
three-dimensional > lent_prepare_study_variants.py -m block -s duineveld.parameter
```

## ALE Relative Reference Frame (RRF)

The ALE relative reference frame cannot use the same boundary conditions because the frame is moving with the droplet. Mesh refinement can be uniform (`-m block` option for `lent_prepare_study_variants.py`) or localy refined (using `cartesianMesh` from the cfMesh project).

Running the rising bubble in the ALE-RRF mode requires the construction of a mesh that only surrounds the near vicinity of the bubble. For this purpose, we start by generating a uniform mesh around the Front, which is then used to reconstruct the Front. The reconstructed Front is then used to create the mesh for the ALE-RRF. A bounding box of the front is computed and a mesh of the bounding box is generated. The Front, reconstructed on the block mesh used for the Eulerian system, is then used to refine the ALE-RRF mesh in the near vicinity of the front. Once the mesh is refined around the interface, new search and distance fields are computed, from which a new front is then reconstructed.

**Note**: The ALE-RRF study requires the change of boundary conditions for the velocity. You will find the required BCs for `U` commented out in `templateCase_risingBubble3D/0.org/U`. 

### Local refinement with cfMesh (cartesianMesh)

At this point, the preparation of the ALE-RRF case study is manual:

1. Uncomment the ALE-RRF boundary conditions for `0.org/U` in `templateCase_risingBubble3D`.
2. Generate the study for the Eulerian mesh first with  
    ```
    three-dimensional > lent_prepare_study_variants.py -m block -s duineveld-cfmesh.parameter
    ```
   the parameters are the same as for the block mesh.
3. In each case of the study, run `cartesianMesh.sh`.
   **ATTENTION**: the utility of `surfaceGenerateBoundingBox` does not resolve intersected region, thus after using this utility, patch generated  from droplet surface should be remove by another utility `removeSurfaceFacets`   
4. In each case of the study, run `rm -rf 0 && cp -r 0.org 0 && lentSetFields`. 
5. In each case of the study, run `lentFoam` to start the simulation.
