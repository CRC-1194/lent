{
    // Solve the momentum equation for the cell centered velocity. 
    solve
    (
        UEqn == lent.cellSurfaceTensionForce(mesh, front)
            + fvc::reconstruct(
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
    );

    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
