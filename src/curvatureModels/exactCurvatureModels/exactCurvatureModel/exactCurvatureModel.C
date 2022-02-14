/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::exactCurvatureModel

SourceFiles
    exactCurvatureModelI.H
    exactCurvatureModel.C
    exactCurvatureModelIO.C

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)
 
Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"


\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H" 

#include "exactCurvatureModel.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(exactCurvatureModel, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void exactCurvatureModel::initializeCellCurvatureBuffer(const fvMesh& mesh) const
{
    if (!cellCurvatureInitialized_)
    {
        cellCurvaturePtr_ = std::shared_ptr<volScalarField>{
            new volScalarField
            {
                IOobject(
                    "cell_curvature", 
                    mesh.time().timeName(), 
                    mesh,
                    IOobject::NO_READ, 
                    IOobject::AUTO_WRITE
                ), 
                mesh, 
                dimensionedScalar(
                    "zero", 
                    dimless/dimLength, 
                    0.0
                )
            }
        };

        cellCurvatureInitialized_ = true;
    }
}

void exactCurvatureModel::updateCurvatureBuffers(const fvMesh& mesh, const triSurfaceFront&) const
{
    // Update cell curvature
    initializeCellCurvatureBuffer(mesh);

    const auto& C = mesh.C(); 

    // Set internal cell centered curvature field
    auto& cellCurvature = *cellCurvaturePtr_;
    forAll(cellCurvature, I)
    {
        cellCurvature[I] = curvatureAtPoint(C[I]);
    }

    // Update face curvature
    auto& faceCurvatureField = *curvatureBuffer(mesh);

    const auto& Cf = mesh.Cf();
    
    forAll(Cf, I)
    {
        faceCurvatureField[I] = curvatureAtPoint(Cf[I]);
    }

    // Set the boundary surface centered curvature field
    auto& surfaceCurvatureBoundaries = faceCurvatureField.boundaryFieldRef(); 
    forAll(surfaceCurvatureBoundaries, I)
    {
        auto& surfaceCurvatureBoundary = surfaceCurvatureBoundaries[I];
        const auto& CfBoundary = Cf.boundaryField()[I];

        forAll(surfaceCurvatureBoundary, J)
        {
            surfaceCurvatureBoundary[J] = curvatureAtPoint(CfBoundary[J]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
exactCurvatureModel::exactCurvatureModel(const dictionary& configDict)
:
    curvatureModel{configDict},
    CurvatureBufferLogic<surfaceScalarField, fvMesh, scalar>{"exact_face_curvature"},
    write_{configDict.lookupOrDefault<Switch>("write", false)}
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
std::shared_ptr<volScalarField> exactCurvatureModel::cellCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& front
) const
{
    if (curvatureRequiresUpdate(mesh))
    {
        updateCurvatureBuffers(mesh, front);
        curvatureUpdated(mesh);
    }

    if (write_)
    {
        cellCurvaturePtr_->write(); 
    }

    return cellCurvaturePtr_;
}

std::shared_ptr<surfaceScalarField> exactCurvatureModel::faceCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& front
) const
{
    if (curvatureRequiresUpdate(mesh))
    {
        updateCurvatureBuffers(mesh, front);
        curvatureUpdated(mesh);
    }

    if (write_)
    {
        curvatureBuffer(mesh)->write(); 
    }

    return curvatureBuffer(mesh);
}

std::shared_ptr<volVectorField> exactCurvatureModel::cellInterfaceNormals(
    const fvMesh& mesh,
    const triSurfaceFront&
) const
{
    std::shared_ptr<volVectorField> cellNormalPtr
    {
        new volVectorField
        {
            IOobject(
                "cell_normals", 
                mesh.time().timeName(), 
                mesh,
                IOobject::NO_READ, 
                IOobject::AUTO_WRITE
            ), 
            mesh, 
            dimensionedVector(
                "zero", 
                dimless, 
                vector(0.0, 0.0, 0.0)
            )
        }
    };

    auto& cellNormals = *cellNormalPtr;

    const auto& C = mesh.C();

    forAll(C, I)
    {
        cellNormals[I] = normalAtPoint(C[I]);
    }

    return cellNormalPtr;
}

vector exactCurvatureModel::normalAtPoint(const point&) const
{
    notImplemented("normalAtPoint");

    return vector{};
}
    

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
