/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::frontTriangleCurvatureModel

SourceFiles
    frontTriangleCurvatureModel.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Curvature model based on the surface tension model described in the
    2012 paper of Tukovic and Jasak

\*---------------------------------------------------------------------------*/

#include "frontTriangleCurvatureModel.H"

#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"

#include "lentCommunication.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontTriangleCurvatureModel, 0);
    addToRunTimeSelectionTable(curvatureModel, frontTriangleCurvatureModel, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void frontTriangleCurvatureModel::computeCurvature(const fvMesh& mesh, const triSurfaceFront& front) const
{
    auto& cn = *curvatureBuffer(front);

    const auto& normalCalculator = normalCalculatorTmp_.ref();
    auto frontVertexNormalsTmp = normalCalculator.vertexNormals(mesh, front);
    auto& n = frontVertexNormalsTmp.ref();

    const auto& faces = front.localFaces();
    const auto& p = front.localPoints();
    const auto& triArea = front.magSf();

    forAll(faces, I)
    {
        const auto& f = faces[I];

        cn[I] = 0.5*(
                        ((p[f[1]] - p[f[0]]) ^ (n[f[1]] + n[f[0]]))
                      + ((p[f[2]] - p[f[1]]) ^ (n[f[2]] + n[f[1]]))
                      + ((p[f[0]] - p[f[2]]) ^ (n[f[0]] + n[f[2]]))
                    ) / triArea[I];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
frontTriangleCurvatureModel::frontTriangleCurvatureModel(const dictionary& configDict)
:
    frontBasedCurvatureModel{configDict},
    normalCalculatorTmp_{
        frontVertexNormalCalculator::New(configDict.subDict("normalCalculator"))
    }
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
