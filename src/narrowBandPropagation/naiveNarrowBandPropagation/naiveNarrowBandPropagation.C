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
    Foam::naiveNarrowBandPropagation

SourceFiles
    naiveNarrowBandPropagation.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)

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

    bool jumpFound = true; 

    while(jumpFound)
    {
        jumpFound = false; 
        forAll (own, faceI)
        {
            if ((signedDistance[own[faceI]] < 0) &&
                (signedDistance[nei[faceI]] == L))
            {
                signedDistance[nei[faceI]] *= -1;
                jumpFound = true; 
            }
            if ((signedDistance[nei[faceI]] < 0) &&
                (signedDistance[own[faceI]] == L))
            {
                signedDistance[own[faceI]] *= -1;
                jumpFound = true; 
            }
        }
    } 

    signedDistance.correctBoundaryConditions(); 
}

void naiveNarrowBandPropagation::ensureNarrowBand(
    pointScalarField& pointSignedDistance,
    scalar L
) const
{
    const pointMesh& pMesh = pointSignedDistance.mesh();

    const labelListList& pointPoints = pMesh().pointPoints();

    //label jumpPoint = -1;
    bool jumpFound = true;

    while(jumpFound)
    {
        jumpFound = false;  

        forAll (pointPoints, pointI)
        {
            const labelList& neighborPoints = pointPoints[pointI];

            forAll(neighborPoints, pointJ)
            {
                if
                (
                    (pointSignedDistance[neighborPoints[pointJ]] == L) &&
                    (pointSignedDistance[pointI] < 0)
                )
                {
                    jumpFound = true; 
                    pointSignedDistance[neighborPoints[pointJ]] *= -1;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
