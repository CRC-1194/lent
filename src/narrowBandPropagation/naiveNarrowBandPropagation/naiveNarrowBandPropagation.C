/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "naiveNarrowBandPropagation.H"
#include "pointMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(naiveNarrowBandPropagation, 0);

    addToRunTimeSelectionTable(narrowBandPropagation, naiveNarrowBandPropagation, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

naiveNarrowBandPropagation::naiveNarrowBandPropagation(const dictionary& configDict) 
:
    narrowBandPropagation(configDict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

naiveNarrowBandPropagation::~naiveNarrowBandPropagation() {}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void naiveNarrowBandPropagation::ensureNarrowBand(
    volScalarField& signedDistance, 
    scalar L
) const
{
    const fvMesh& mesh = signedDistance.mesh(); 
    const labelList& own = mesh.owner(); 
    const labelList& nei = mesh.neighbour(); 

    label jumpFace = -1; 

    do 
    {
        jumpFace = -1; 
        forAll (own, faceI)
        {
            if ((signedDistance[own[faceI]] < 0) && 
                (signedDistance[nei[faceI]] == L))
            {
                jumpFace = faceI; 
                signedDistance[nei[faceI]] *= -1; 
            }
            if ((signedDistance[nei[faceI]] < 0) &&
                (signedDistance[own[faceI]] == L))
            {
                jumpFace = faceI; 
                signedDistance[own[faceI]] *= -1;
            }
        }
    } while (jumpFace >= 0);

    signedDistance.boundaryField().evaluate(); 
}

void naiveNarrowBandPropagation::ensureNarrowBand(
    pointScalarField& pointSignedDistance, 
    scalar L
) const
{
    const pointMesh& pMesh = pointSignedDistance.mesh(); 

    const labelListList& pointPoints = pMesh().pointPoints(); 

    label jumpPoint = -1; 

    do 
    {
        jumpPoint = -1; 

        forAll (pointPoints, pointI)
        {
            const labelList& pointIpoints = pointPoints[pointI];  
            
            forAll(pointIpoints, pointJ)
            {
                if 
                (
                    (pointSignedDistance[pointIpoints[pointJ]] == L) && 
                    (pointSignedDistance[pointI] < 0)
                ) 
                {
                    jumpPoint = pointIpoints[pointJ]; 
                    pointSignedDistance[pointIpoints[pointJ]] *= -1; 
                }
            }
        }

    } while (jumpPoint >= 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
