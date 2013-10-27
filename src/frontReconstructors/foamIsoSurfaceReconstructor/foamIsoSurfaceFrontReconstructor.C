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

#include "foamIsoSurfaceFrontReconstructor.H"
#include "addToRunTimeSelectionTable.H"
#include "isoSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(foamIsoSurfaceFrontReconstructor, 0); 
    addToRunTimeSelectionTable(foamIsoSurfaceFrontReconstructor, foamIsoSurfaceFrontReconstructor, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

foamIsoSurfaceFrontReconstructor::foamIsoSurfaceFrontReconstructor(
   const dictionary& configDict
)
:
    frontReconstructor(configDict), 
    mergeTolerance_(readScalar(configDict.lookup("mergeTolerance"))), 
    regularize_(configDict.lookup("regularization"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

foamIsoSurfaceFrontReconstructor::~foamIsoSurfaceFrontReconstructor()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

labelList foamIsoSurfaceFrontReconstructor::reconstructFront(
    triSurfaceFront& front, 
    const volScalarField& signedDistance,
    const pointScalarField& pointSignedDistance
) const
{
    isoSurface iso (
        signedDistance, 
        pointSignedDistance, 
        0, 
        regularize_,
        mergeTolerance_
    );

    front = iso; 

    return iso.meshCells(); 
} 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
