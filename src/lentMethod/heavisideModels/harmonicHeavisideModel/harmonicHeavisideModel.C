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

#include "harmonicHeavisideModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 
    
    defineTypeNameAndDebug(harmonicHeavisideModel, 0); 
    addToRunTimeSelectionTable(heavisideModel, harmonicHeavisideModel, Empty);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

harmonicHeavisideModel::harmonicHeavisideModel()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

harmonicHeavisideModel::~harmonicHeavisideModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void harmonicHeavisideModel::calcHeavisideField(
    volScalarField& heaviside, 
    const volScalarField& signedDistance, 
    const volScalarField& narrowBandWidthSqr 
) const
{
    scalar pi = constant::mathematical::pi; 

    forAll (heaviside, cellI)
    {
        scalar narrowBandWidth = sqrt(narrowBandWidthSqr[cellI]);

        if (mag(signedDistance[cellI]) < narrowBandWidth)
        {
            heaviside[cellI] = 0.5 * (
                1 + signedDistance[cellI] / narrowBandWidth + 1/pi * 
                sin((pi * signedDistance[cellI]) / narrowBandWidth)
            );
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

