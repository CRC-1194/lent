{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    fvOptions.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            lent.faceSurfaceTensionForce(mesh, frontMesh)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()

    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();
            // Update the momentum flux explictly. 
            p_rgh.relax();

            // FIXME: Test the volume conservation with the two ways of
            // correcting the velocity. TM
            
            // First order accurate explicit update of the velocity field.
            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            
            // Assemble the matrix anew. TM.  
            // Update the momentum flux using the volume conservative flux.
            rhoPhi == fvc::interpolate(rho) * phi; 
            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(rhoPhi, U)
              + turbulence->divDevRhoReff(rho, U)
             ==
                fvOptions(rho, U)
            );
            UEqn.relax();
            fvOptions.constrain(UEqn);

            // Solve the momentum equaiton for the cell centered velocity. 
            solve
            (
                UEqn == lent.cellSurfaceTensionForce(mesh, frontMesh)
                    + fvc::reconstruct(
                        (
                          - ghf*fvc::snGrad(rho)
                          - fvc::snGrad(p_rgh)
                        ) * mesh.magSf()
                    )
                    - fvc::grad(p_rgh)  // Add the final pressure to the equation.
            );

            U.correctBoundaryConditions();
            fvOptions.correct(U);
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
