    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      // see method divDevRhoReff(...) in "linearViscousStress.C"
      // for OpenFOAM's discretization of the diffusion term 
      + turbulence->divDevRhoReff(rho, U)
      - lent.surfaceTensionImplicitPart(U, markerField, front)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (lentSC.momentumPredictor())
    {
        solve
        (
            UEqn
         ==

            lent.cellSurfaceTensionForce(mesh, front)

            +

            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
