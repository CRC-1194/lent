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

#include "heavisideModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(heavisideModel, 0); 
    defineRunTimeSelectionTable(heavisideModel, Empty);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

heavisideModel::heavisideModel()
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<heavisideModel>
heavisideModel::New(const word& name)
{
    if (debug)
    {
        Info<< "Selecting heavisideModel" << name << endl;
    }

    // Find the constructor pointer for the model in the constructor table.
    EmptyConstructorTable::iterator cstrIter =
        EmptyConstructorTablePtr_->find(name);

    // If the constructor pointer is not found in the table.
    if (cstrIter == EmptyConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "heavisideModel::New(const dictionary&)"
        )   << "Unknown heavisideModel type "
            << name << nl << nl
            << "Valid heavisideModels are : " << endl
            << EmptyConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // Construct the model and return the autoPtr to the object. 
    return autoPtr<heavisideModel>
            (cstrIter()(name));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heavisideModel::~heavisideModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool heavisideModel::distanceWithinNarrowBand(
    scalar distance, 
    scalar narrowBandWidth
) const
{
    return mag(distance) < narrowBandWidth;  
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void heavisideModel::operator=(const heavisideModel& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("heavisideModel::operator=(const heavisideModel&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
