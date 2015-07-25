scalar computeU(vector v);

scalar computeV(vector v);

scalar ellipsoidDeviation(scalar u, scalar v, vector axis, vector& p,
                          vector center);

vector ellipsoidNormal(scalar u, scalar v, vector axis);

scalar ellipsoidCurvature(scalar u, scalar v, vector axis);

void checkEllipsoidCurvature(const triSurfacePointVectorField& cn,
                             const triSurface& front,
                             const triSurfacePointVectorField& uv,
                             vector center,
                             vector axis, std::fstream& errorFile);

void checkEllipsoidNormal(const triSurfacePointVectorField& cn,
                          const triSurface& front,
                          const triSurfacePointVectorField& uv,
                          vector center,
                          vector axis, std::fstream& errorFile);

void ellipsoidDeviation(const triSurface& front,
                        const triSurfacePointVectorField& uv,
                        vector center, vector axis,
                        std::fstream& errorFile);

void correctFront(triSurface& front,
                  triSurfacePointVectorField& uv,
                  vector center, vector axis);
