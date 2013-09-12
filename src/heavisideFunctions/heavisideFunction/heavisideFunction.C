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

#include "heavisideFunction.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
    defineTypeNameAndDebug(heavisideFunction, 0); 
    defineRunTimeSelectionTable(heavisideFunction, Empty);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heavisideFunction::heavisideFunction()
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::heavisideFunction>
Foam::heavisideFunction::New(const word& name)
{
    if (debug)
    {
        Info<< "Selecting heavisideFunction" << name << endl;
    }

    // Find the constructor pointer for the model in the constructor table.
    EmptyConstructorTable::iterator cstrIter =
        EmptyConstructorTablePtr_->find(name);

    // If the constructor pointer is not found in the table.
    if (cstrIter == EmptyConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "heavisideFunction::New(const dictionary&)"
        )   << "Unknown heavisideFunction type "
            << name << nl << nl
            << "Valid heavisideFunctions are : " << endl
            << EmptyConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // Construct the model and return the autoPtr to the object. 
    return autoPtr<heavisideFunction>
            (cstrIter()(name));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heavisideFunction::~heavisideFunction()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::heavisideFunction::distanceWithinNarrowBand(
    scalar distance, 
    scalar narrowBandWidth
) const
{
    return mag(distance) < narrowBandWidth;  
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::heavisideFunction::operator=(const heavisideFunction& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::heavisideFunction::operator=(const Foam::heavisideFunction&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// ************************************************************************* //
