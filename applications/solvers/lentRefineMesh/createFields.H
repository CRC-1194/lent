
    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh(
        IOobject(
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

    Info << "Reading the cell-centered distance field signedDistance\n" << endl;
    volScalarField signedDistance(
        IOobject(
            "signedDistance",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField searchDistanceSqr 
    (
        IOobject
        (
            "searchDistanceSqr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            // TODO: uncomment, used only for debugging
            //IOobject::NO_WRITE,
            IOobject::AUTO_WRITE
        ),
        mesh, 
        dimensionedScalar(
            "zero",
            dimLength,
            0
        ), 
        "zeroGradient"
    );

    pointMesh pmesh(mesh); 

    // Initialize the points-to-elements distance. 
    pointScalarField pointSignedDistance(
        IOobject(
            "pointSignedDistance", 
            runTime.timeName(), 
            mesh, 
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pmesh
    );

    pointScalarField pointSearchDistanceSqr
    (
        IOobject
        (
            "pointSearchDistanceSqr", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        pmesh,
        dimensionedScalar(
            "zero",
            dimLength,
            0
        ),
        "fixedValue"
    );

    Info << "Reading the cell-centered Heaviside function field H\n" << endl;
    volScalarField heaviside(
        IOobject(
            "heaviside",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    incompressibleTwoPhaseMixture twoPhaseProperties(U, phi, "heaviside");

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho(
        IOobject(
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        heaviside*rho1 + (scalar(1) - heaviside)*rho2,
        heaviside.boundaryField().types()
    );

    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi(
        IOobject(
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    // Construct interface from alpha1 distribution
    interfaceProperties interface(heaviside, U, twoPhaseProperties);
