{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
        // TODO: find out what this correction does. If it yields a non-zero
        // contribution, fix the interpolation of rho as in the new calculation
        // of "rhoPhi"
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            lent.faceSurfaceTensionForce(mesh, front)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (lentSC.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(lentSC.finalInnerIter())));

        if (lentSC.finalNonOrthogonalIter())
        {
            // Update of volumetric fluxes
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            // First order accurate explicit update of the velocity field.
            if (lentSC.explicitVelocityUpdate())
            {
                #include "U_reconstruct.H"
            }
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
