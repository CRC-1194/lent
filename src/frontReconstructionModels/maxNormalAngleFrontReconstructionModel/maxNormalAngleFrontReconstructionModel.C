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
    Foam::maxNormalAngleFrontReconstructionModel

SourceFiles
    maxNormalAngleFrontReconstructionModel.C

Author
    Tomislav Maric (maric@csi.tu-darmstadt.de)

Description

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


#include "maxNormalAngleFrontReconstructionModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(maxNormalAngleFrontReconstructionModel, 0);
    addToRunTimeSelectionTable(frontReconstructionModel, maxNormalAngleFrontReconstructionModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

maxNormalAngleFrontReconstructionModel::maxNormalAngleFrontReconstructionModel(const dictionary& configDict)
:
    frontReconstructionModel(configDict),
    maxAngle_(configDict.get<scalar>("value") * M_PI / 180.0),
    minAngleCos_(Foam::cos(maxAngle_)),
    previouslyReconstructed_(false)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool maxNormalAngleFrontReconstructionModel::reconstructionRequired(
    const triSurfaceFront& front,
    const volScalarField& signedDistance
) const
{
    if (signedDistance.time().timeIndex() <= 1)
        return true; 

    const auto& edges = front.edges(); 
    const auto& allEdgeFaces = front.edgeFaces(); 
    const auto& faceNormals = front.faceNormals(); 

    forAll(edges, I)
    {
        const auto& edgeFaces = allEdgeFaces[I];
        const auto& n0 = faceNormals[edgeFaces[0]]; 

        for(label J = 1; J < edgeFaces.size(); ++J)
        {
            const auto& n = faceNormals[edgeFaces[J]]; 

            if (((n0 & n) < minAngleCos_) /*&& (! previouslyReconstructed_)*/)
            {
                previouslyReconstructed_ = true; 
                return true; 
            }
        }
    }

    previouslyReconstructed_ = false; 

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

