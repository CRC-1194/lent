    // The fields defined in this file are required for any test class
    // which is derived from "lentSubalgorithmTest".
    // So include it the corresponding test application.

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
