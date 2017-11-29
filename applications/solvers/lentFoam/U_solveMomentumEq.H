{
    // Do nit update the U sytem of equations here as it will break consistency
    // with the p equation
    //rhoPhi == fvc::interpolate(rho) * phi;
    //#include "computeRhoPhi.H"
    /*
    UEqn =
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);
    */

    // Solve the momentum equation for the cell centered velocity. 
    solve
    (
        UEqn == lent.cellSurfaceTensionForce(mesh, front)
                    - fvc::grad(p_rgh)
                    - gh*fvc::grad(rho)
            //+ fvc::reconstruct(
            //    (
            //      - ghf*fvc::snGrad(rho)
            //      - fvc::snGrad(p_rgh)
            //    ) * mesh.magSf()
            //)
    );

    U.correctBoundaryConditions();
    fvOptions.correct(U);
}