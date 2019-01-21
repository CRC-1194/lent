    // Initialize the test fields.
    // - Ellipsoid
    auto veFieldTmp = makeField<volScalarField>("veField", mesh);  
    auto& veField = veFieldTmp.ref(); 
    pointMesh pMesh(mesh); 
    pointScalarField peField
    (
        IOobject
        (
            "peField", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        pMesh, 
        dimensionedScalar("pField", dimless, 0)
    );
    // - Sphere 
    auto vsFieldTmp = makeField<volScalarField>("vsField", mesh);  
    auto& vsField = vsFieldTmp.ref(); 
    pointScalarField psField
    (
        IOobject
        (
            "psField", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        pMesh, 
        dimensionedScalar("pField", dimless, 0)
    );

    // Initialie error fields.
    // - Ellipsoid
    auto linfInit = dimensionedScalar("Linf", dimLength, -VGREAT);
    auto LinfSphereBccTmp = makeField<volScalarField>("LinfSphereBcc", mesh, linfInit); 
    auto& LinfSphereBcc = LinfSphereBccTmp.ref();
    auto LinfSphereBccFvmTmp = makeField<volScalarField>("LinfSphereBccFvm", mesh, linfInit); 
    auto& LinfSphereBccFvm = LinfSphereBccFvmTmp.ref();
    // - Sphere 
    auto LinfEllipsoidBccTmp = makeField<volScalarField>("LinfEllipsoidBcc", mesh, linfInit); 
    auto& LinfEllipsoidBcc = LinfEllipsoidBccTmp.ref();
    auto LinfEllipsoidBccFvmTmp = makeField<volScalarField>("LinfEllipsoidBccFvm", mesh, linfInit); 
    auto& LinfEllipsoidBccFvm = LinfEllipsoidBccFvmTmp.ref();

    // Set ellipsoid and sphere properties for convergence tests.
    std::vector<double> aAxes = {0.046875, 0.09375, 0.1875, 0.375, 0.75}; 
    std::vector<double> bAxes(aAxes);
    std::vector<double> cAxes(aAxes);
    for (auto& axis : bAxes)
        axis *= 2.0 / 3.0; 
    for (auto& axis : cAxes)
        axis *= 1.0 / 3.0; 
    std::vector<double> radii = {0.03125,0.0625,0.125,0.25,0.5};
    std::vector<int> cellsPerLength = {4,8,16,32,64}; 
    assert(radii.size() == aAxes.size() && bAxes.size() == cAxes.size()); 

    using pointVectorVector = std::vector<std::vector<point>>;
    const label nPointsPerCell = 100; 
    pointVectorVector randPointsInCells
    (
        mesh.nCells(), 
        std::vector<point>(nPointsPerCell, point(VGREAT, VGREAT, VGREAT))
    );
    Info << "Generating " << nPointsPerCell << " random sampling points per cell ..."; 
    genRandomPointsInCells(randPointsInCells, mesh, nPointsPerCell);
    Info << "Done." << endl;
