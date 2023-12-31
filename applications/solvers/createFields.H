    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading the cell distance field signedDistance\n" << endl;
    volScalarField signedDistance
    (
        IOobject
        (
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
            IOobject::NO_WRITE
        ),
        mesh, 
        dimensionedScalar(
            "zero",
            dimLength,
            0
        )
    );

    pointMesh pmesh(mesh); 

    // Initialize the points-to-elements distance. 
    Info << "Reading the point distance field signedDistance\n" << endl;
    pointScalarField pointSignedDistance 
    (
        IOobject
        (
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
            IOobject::NO_WRITE
        ),
        pmesh,
        dimensionedScalar(
            "zero",
            dimLength,
            0
        )
    );

    #include "createPhi.H"

    // FIXME: check if this can be avoided. The compilation 
    // of the incompressibleTwoPhaseMixture holds relative 
    // paths from within the folder with `Make`. It makes 
    // it impossible to be compiled into any other code than
    // the mixture library. Hence I can't lookup mixture from
    // the registry - it is using multiple inheritance from 
    // transportModel, IOdictionary and twoPhaseMixture. TM.
    IOdictionary transportProperties 
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            runTime, 
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ); 

    incompressibleTwoPhaseMixture mixture(U, phi);
    volScalarField& markerField(mixture.alpha1());

    
    // FIXME: Clean up: remove mixture and its trainwrecks. TM.
    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar nu1 = mixture.nuModel1().viscosityProperties().get<dimensionedScalar>("nu");

    const dimensionedScalar& rho2 = mixture.rho2();
    const dimensionedScalar nu2 = mixture.nuModel2().viscosityProperties().get<dimensionedScalar>("nu");

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        markerField*rho1 + (scalar(1) - markerField)*rho2,
        markerField.boundaryField().types()
    );

    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );

    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, mixture)
    );

    // Construct interface from markerField distribution
    //interfaceProperties interface(markerField, U, mixture);

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

    volScalarField p
    (
        IOobject
        (
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
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("lentSC"),
        pRefCell,
        pRefValue
    );

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

