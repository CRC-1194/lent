#include "fvCFD.H"
#include "triSurface.H"
#include "triSurfaceFields.H"

#include "frontCurvatureModel.H"

// Constructors
frontCurvatureModel::frontCurvatureModel(const fileName& frontName)
{
    triSurface front_(frontName);
    curvatureNormals = dimensionedVector("zero",
                            dimless/dimLength,
                            vector(0,0,0)
                           );
}

// Destructors

// Private member functions

// Public member functions
const triSurfacePointVectorfield& frontCurvatureModel::curvatureNormals()
{
    return curvatureNormals_;
}

void frontCurvatureModel::curvature(triSurfacePointScalarField& curvature)
{
    curvature = mag(curvatureNormals);
}

void frontCurvatureModel::normal(triSurfacePointVectorField& normal)
{
    normal = curvatureNormals / mag(curvatureNormals);
}

void frontCurvatureModel::update()
{
    calcCurvatureNormals();
}
