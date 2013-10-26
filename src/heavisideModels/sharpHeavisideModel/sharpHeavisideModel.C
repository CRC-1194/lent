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

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#include "sharpHeavisideModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 
    
    defineTypeNameAndDebug(sharpHeavisideModel, 0); 
    addToRunTimeSelectionTable(heavisideModel, sharpHeavisideModel, Empty);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sharpHeavisideModel::sharpHeavisideModel()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sharpHeavisideModel::~sharpHeavisideModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void sharpHeavisideModel::calcHeavisideField(
    volScalarField& heaviside, 
    const volScalarField& signedDistance, 
    const volScalarField& searchDistanceSqr 
) const
{
    Info << "SHARP HEAVISIDE FIELD" << endl;

    forAll (heaviside, cellI)
    {
        scalar searchDistance = sqrt(searchDistanceSqr[cellI]);

        if (mag(signedDistance[cellI]) < searchDistance)
        {
            heaviside[cellI] = 0.5;
        } 
        else
        {
            if (signedDistance[cellI] > 0)
            {
                heaviside[cellI] = 1; 
            }
            if (signedDistance[cellI] < 0)
            {
                heaviside[cellI] = 0;
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

