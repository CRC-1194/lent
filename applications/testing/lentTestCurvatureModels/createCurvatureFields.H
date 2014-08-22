Info << "Reading exact face centered curvature..." << endl;
surfaceScalarField faceCurvatureExact
(
    IOobject
    (
        "faceCurvatureExact",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh 
);

Info << "Reading exact cell centered curvature..." << endl;
volScalarField cellCurvatureExact
(
    IOobject
    (
        "cellCurvature",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh 
);
