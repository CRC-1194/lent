
    Info<< "Reading field U\n" << endl;
    volVectorField U(
        IOobject(
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading the cell-centered MarkerField function field H\n" << endl;
    volScalarField markerField(
        IOobject(
            "markerField",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;
    incompressibleTwoPhaseMixture twoPhaseProperties(U, phi, "markerField");

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

    // Construct interface from alpha1 distribution
    interfaceProperties interface(markerField, U, twoPhaseProperties);

