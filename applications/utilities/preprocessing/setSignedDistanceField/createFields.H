    Info << "Reading the cells-to-elements distance field signedDistance\n"; 
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
    Info << "Done." << endl;

    Info << "Initializing points-to-elements distance field pointSignedDistance\n"; 
    pointScalarField pointSignedDistance 
    (
        IOobject
        (
            "pointSignedDistance", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        pointMesh(mesh), 
        GREAT,
        "zeroGradient"
    );
    Info << "Done." << endl;
