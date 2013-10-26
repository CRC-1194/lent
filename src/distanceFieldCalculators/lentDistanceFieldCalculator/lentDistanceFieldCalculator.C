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

#include "lentDistanceFieldCalculator.H"

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(lentDistanceFieldCalculator, 0);
    defineRunTimeSelectionTable(lentDistanceFieldCalculator, Dictionary);

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentDistanceFieldCalculator::lentDistanceFieldCalculator(const dictionary& config)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


tmp<lentDistanceFieldCalculator>
lentDistanceFieldCalculator::New(
   const word& name,
   const dictionary& configDict
)
{
    if (debug)
    {
        Info<< "Selecting lentDistanceFieldCalculator" << name << endl;
    }

    // Find the constructor pointer for the model in the constructor table.
    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    // If the constructor pointer is not found in the table.
    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "lentDistanceFieldCalculator::New(const word& name)"
        )   << "Unknown lentDistanceFieldCalculator type "
            << name << nl << nl
            << "Valid lentDistanceFieldCalculators are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // Construct the model and return the autoPtr to the object. 
    return tmp<lentDistanceFieldCalculator> (cstrIter()(configDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lentDistanceFieldCalculator::~lentDistanceFieldCalculator()
{}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
