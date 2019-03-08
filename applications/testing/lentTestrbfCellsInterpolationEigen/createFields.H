    // Initialize the surface test fields.
    surfaceTestFields ellFields(mesh, "Ellipsoid"); 
    surfaceTestFields sphFields(mesh, "Sphere"); 

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
