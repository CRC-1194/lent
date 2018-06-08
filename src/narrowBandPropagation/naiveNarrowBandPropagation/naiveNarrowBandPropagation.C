/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
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

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

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
