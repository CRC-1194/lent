/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 AUTHOR,AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "frontBasedCurvatureModel.H"
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
