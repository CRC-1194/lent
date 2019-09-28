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

#include "averagingCurvatureModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcAverage.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(averagingCurvatureModel, 0);
    addToRunTimeSelectionTable(curvatureModel, averagingCurvatureModel, Dictionary);
    
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void averagingCurvatureModel::computeCurvature(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const auto& curvatureInputField = inputField(mesh);

    // Compute the curvature using an averaged field.
    volScalarField inputFieldSmooth("smoothInputField", curvatureInputField); 

    for (label I = 0; I < averagingIterations_; ++I)
    {
        inputFieldSmooth == fvc::average(inputFieldSmooth);
    }

    auto cellCurvatureTmp = levelSetCurvature(inputFieldSmooth);
    auto curvatureBufferPtr = curvatureBuffer(mesh);

    // TODO: performance optimization: avoid copying of the curvature field (TT)
    *curvatureBufferPtr = cellCurvatureTmp.ref();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
averagingCurvatureModel::averagingCurvatureModel(const dictionary& configDict)
:
    divergenceBasedCurvatureModel{configDict},
    averagingIterations_{readLabel(configDict.lookup("averagingIterations"))}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
