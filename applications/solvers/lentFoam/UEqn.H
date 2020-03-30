    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U) 
      - fvm::div(rhof*(Ub & mesh.Sf()), U)
      + MRF.DDt(rho, U)
      // see method divDevRhoReff(...) in "linearViscousStress.C"
      // for OpenFOAM's discretization of the diffusion term 
      // TODO: Remove, old code
      //+ turbulence->divDevRhoReff(rho, U)
      - fvm::laplacian(muf, U)
      - fvc::div(muf * (fvc::interpolate(dev2(T(fvc::grad(U)))) & mesh.Sf()))
      // Semi-implicit handling of surface tension forces 
      // Hysing, S. (2006). A new implicit surface tension implementation for
      // interfacial flows. International Journal for Numerical Methods in Fluids,
      // 51(6), 659¿672. https://doi.org/10.1002/fld.1147
      // - lent.surfaceTensionImplicitPart(U, markerField, front)
     ==
        fvOptions(rho, U)
    );
    

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (lentSC.momentumPredictor())
    {
        fvVectorMatrix UEqnWithRhs{
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
        };
        lentSC.solveForVelocity(UEqnWithRhs);

        fvOptions.correct(U);
    }

