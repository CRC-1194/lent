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

#include "normalConsistency.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(normalConsistency, 0); 
    defineRunTimeSelectionTable(normalConsistency, Dictionary);
    addToRunTimeSelectionTable(normalConsistency, normalConsistency, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

normalConsistency::normalConsistency(const dictionary& configDict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void normalConsistency::makeFrontNormalsConsistent(
    triSurfaceFront& front, 
    const volScalarField& signedDistance
) const
{
    // Gradient based normal consistency algorithm.
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<normalConsistency>
normalConsistency::New(const dictionary& configDict)
{
    const word name = configDict.lookup("type"); 

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "normalConsistency::New(const word& name)"
        )   << "Unknown normalConsistency type "
            << name << nl << nl
            << "Valid normalConsistencys are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<normalConsistency> (cstrIter()(configDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

normalConsistency::~normalConsistency()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
