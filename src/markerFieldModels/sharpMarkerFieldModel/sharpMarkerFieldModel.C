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

Class
    Foam::sharpMarkerFieldModel

Description
    Abstract base class for the markerField function calculation from a signed 
    distance field.

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "sharpMarkerFieldModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 
    
    defineTypeNameAndDebug(sharpMarkerFieldModel, 0); 
    addToRunTimeSelectionTable(markerFieldModel, sharpMarkerFieldModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sharpMarkerFieldModel::sharpMarkerFieldModel(const dictionary& configDict)
:
    markerFieldModel(configDict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sharpMarkerFieldModel::~sharpMarkerFieldModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void sharpMarkerFieldModel::calcMarkerField(
    volScalarField& markerField, 
    const volScalarField& signedDistance,
    const volScalarField& searchDistanceSqr
) const
{
    forAll (markerField, cellI)
    {
        scalar searchDistance = sqrt(searchDistanceSqr[cellI]);

        if (mag(signedDistance[cellI]) < searchDistance)
        {
            markerField[cellI] = 0.5;
        } 
        else
        {
            if (signedDistance[cellI] > 0)
            {
                markerField[cellI] = 1; 
            }
            if (signedDistance[cellI] < 0)
            {
                markerField[cellI] = 0;
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

