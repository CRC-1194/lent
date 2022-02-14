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
    Foam::frontBasedCurvatureModel

SourceFiles
    frontBasedCurvatureModelI.H
    frontBasedCurvatureModel.C
    frontBasedCurvatureModelIO.C

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

#include "frontBasedCurvatureModel.H"
#include "lentCommunication.H"
#include "triSurfaceFrontFields.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontBasedCurvatureModel, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void frontBasedCurvatureModel::updateCurvatureBuffers(const fvMesh& mesh, const triSurfaceFront& front) const
{
    if (curvatureRequiresUpdate(mesh))
    {
        resizeBufferField(front.UList<labelledTri>::size(), front);
        setBufferToZero(front);
        computeCurvature(mesh, front);
        frontToMeshTmp_->transferCurvature(*curvatureBuffer(front), front, mesh); 
        curvatureUpdated(mesh);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
frontBasedCurvatureModel::frontBasedCurvatureModel(const dictionary& configDict)
:
    curvatureModel{configDict},
    CurvatureBufferLogic<triSurfaceFrontVectorField, triSurfaceFront, vector>{"front_curvature_normal"},
    frontToMeshTmp_{frontToMeshTransferModel::New(configDict.subDict("frontToMeshTransfer"))}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
std::shared_ptr<volScalarField> frontBasedCurvatureModel::cellCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& front
) const
{
    updateCurvatureBuffers(mesh, front);
    return frontToMeshTmp_->cellCurvature(mesh);
}

std::shared_ptr<surfaceScalarField> frontBasedCurvatureModel::faceCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& front
) const 
{
    updateCurvatureBuffers(mesh, front);
    return frontToMeshTmp_->faceCurvature(mesh);
}

std::shared_ptr<triSurfaceFrontVectorField> frontBasedCurvatureModel::frontCurvatureNormal(
    const fvMesh& mesh,
    const triSurfaceFront& front
) const
{
    updateCurvatureBuffers(mesh, front);
    return curvatureBuffer(front);
}

std::shared_ptr<volVectorField> frontBasedCurvatureModel::cellInterfaceNormals(
    const fvMesh& mesh,
    const triSurfaceFront& front
) const
{
    // TODO: this is just a preliminary approach to cell interface normal
    // computation (TT)
    std::shared_ptr<volVectorField> cellNormalPtr{
        new volVectorField{
            IOobject
            (
                "cell_interface_normals",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector{
                "zero",
                dimless,
                vector{0, 0, 0}
            }
        }
    };

    auto& cellNormals = *cellNormalPtr;

    // TODO: tranfer logic is taken from "triangleInCellTranferModel" and
    // duplicated here. Remove duplication when refactoring (TT)
    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 

    const auto& trianglesInCell = communication.interfaceCellToTriangles();
    const auto& faceNormal = front.Sf();

    // Area weighted averaging of triangle normals
    for (const auto& cellTrianglesMap : trianglesInCell)
    {
        const auto& cellLabel = cellTrianglesMap.first;
        const auto& triangleLabels = cellTrianglesMap.second;

        for (const auto& tl : triangleLabels)
        {
            cellNormals[cellLabel] += faceNormal[tl];
        }

        cellNormals[cellLabel] /= mag(cellNormals[cellLabel]) + SMALL;
    }

    // Narrow band propagation
    const auto& cellToTriangle = communication.cellsTriangleNearest();
    const auto& triangleToCell = communication.triangleToCell();

    forAll(cellToTriangle, I)
    {
        const auto& hitObject = cellToTriangle[I];
        
        if (hitObject.hit())
        {
            cellNormals[I] = cellNormals[triangleToCell[hitObject.index()]];
        }
    }

    return cellNormalPtr;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
