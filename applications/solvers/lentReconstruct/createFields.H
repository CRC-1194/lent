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
            IOobject::AUTO_WRITE
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

    Info << "Reading the cell-centered MarkerField field markerField\n" << endl;
    volScalarField markerField 
    (
        IOobject
        (
            "markerField",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh
    );
