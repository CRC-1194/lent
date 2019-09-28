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

#include "addToRunTimeSelectionTable.H" 
#include "fvcGrad.H" 
#include "fvcDiv.H"
#include "surfaceInterpolate.H" 

#include "divergenceBasedCurvatureModel.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(divergenceBasedCurvatureModel, 0);
    addToRunTimeSelectionTable(curvatureModel, divergenceBasedCurvatureModel, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void divergenceBasedCurvatureModel::computeCurvature(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const auto& curvatureInputField = inputField(mesh);

    auto cellCurvatureTmp = levelSetCurvature(curvatureInputField);
    auto curvatureBufferPtr = curvatureBuffer(mesh);

    // TODO: performance optimization: avoid copying of the curvature field (TT)
    *curvatureBufferPtr = cellCurvatureTmp.ref();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
tmp<volScalarField> divergenceBasedCurvatureModel::levelSetCurvature(const volScalarField& levelSetField) const
{
    const surfaceVectorField& Sf = levelSetField.mesh().Sf();

    //Cell gradient of levelSetField
    const volVectorField curvGrad(fvc::grad(levelSetField, "curvatureGradient"));

    // Interpolated face-gradient of levelSetField 
    surfaceVectorField curvGradF(fvc::interpolate(curvGrad));

    // Hardcoded stabilization of the gradient to avoid floating point
    // exception.
    dimensionedScalar deltaN
    (
        "deltaN",
        levelSetField.dimensions() / dimLength, 
        SMALL 
    );

    // Face unit interface normal
    surfaceVectorField curvGradFhat(curvGradF /(mag(curvGradF) + deltaN));

    return -fvc::div(curvGradFhat & Sf); 
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
divergenceBasedCurvatureModel::divergenceBasedCurvatureModel(const dictionary& configDict)
:
    curvatureModel{configDict},
    CurvatureBufferLogic<volScalarField, fvMesh, scalar>{"cell_curvature"},
    inputFieldName_{configDict.lookup("curvatureField")}
{
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
std::shared_ptr<volScalarField> divergenceBasedCurvatureModel::cellCurvature(
        const fvMesh& mesh, 
        const triSurfaceFront& front
    ) const
{
    if (curvatureRequiresUpdate(mesh))
    {
        computeCurvature(mesh, front);
        curvatureUpdated(mesh);
    }

    return curvatureBuffer(mesh);
}

std::shared_ptr<surfaceScalarField> divergenceBasedCurvatureModel::faceCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& front
) const
{
    auto& cellCurvatureField = *cellCurvature(mesh, front);

    auto faceCurvatureTmp = fvc::interpolate(cellCurvatureField);

    return std::shared_ptr<surfaceScalarField>{faceCurvatureTmp.ptr()};
}

std::shared_ptr<volVectorField> divergenceBasedCurvatureModel::cellInterfaceNormals(
    const fvMesh& mesh,
    const triSurfaceFront& front
) const
{
    const auto& curvatureInputField = inputField(mesh);

    auto interfaceNormalsTmp = fvc::grad(curvatureInputField)/(mag(fvc::grad(curvatureInputField)) + SMALL);

    return std::shared_ptr<volVectorField>{interfaceNormalsTmp.ptr()};
}
    
const word& divergenceBasedCurvatureModel::inputFieldName() const
{
    return inputFieldName_;
}

const volScalarField& divergenceBasedCurvatureModel::inputField(const fvMesh& mesh) const
{
    return mesh.lookupObject<volScalarField>(inputFieldName()); 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
