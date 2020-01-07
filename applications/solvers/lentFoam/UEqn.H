    // The momentum flux is computed from MULES as  
    // rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2; 
    // However, LENT has no ability to compute the volumetric phase flux. 
    // TODO: examine the impact of the momentum flux computation and devise
    // more accurate approach if required (TT)
    rhoPhi.oldTime() == rhoPhi;
    rhoPhi == rhof * phi;

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvc::div(rhoPhi, U) 
      - fvm::div(rhof*(Ub & mesh.Sf()), U)
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

