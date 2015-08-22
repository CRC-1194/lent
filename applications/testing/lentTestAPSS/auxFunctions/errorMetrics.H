void checkCurvature(const triSurfacePointVectorField& cn, const triSurface& front,
                    scalar radius, std::fstream& errorFile);

void checkNormal(const triSurfacePointVectorField& cn, const triSurface& front,
                 vector center, std::fstream& errorFile);

void sphereDeviation(const triSurface& front, scalar radius, vector center,
                     std::fstream& errorFile);

void meshQuality(const triSurface& front, std::fstream& file);
