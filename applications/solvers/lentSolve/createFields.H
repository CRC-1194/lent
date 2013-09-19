
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

    // Initialize the points-to-elements distance. 
    pointScalarField pointSignedDistance(
        IOobject(
            "pointSignedDistance", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        pointMesh(mesh), 
        GREAT,
        "zeroGradient"
    );

    Info << "Reading the cell-centered Heaviside field\n" << endl;
    volScalarField H(
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
    incompressibleTwoPhaseMixture twoPhaseProperties(U, phi, "H");

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
        H*rho1 + (scalar(1) - H)*rho2,
        H.boundaryField().types()
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
    interfaceProperties interface(H, U, twoPhaseProperties);


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence(
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    #include "readGravitationalAcceleration.H"

    /*
    dimensionedVector g0(g);

    // Read the data file and initialise the interpolation table
    interpolationTable<vector> timeSeriesAcceleration
    (
        runTime.path()/runTime.caseConstant()/"acceleration.dat"
    );
    */

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p(
        IOobject(
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference()) {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
