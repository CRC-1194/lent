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
    Foam::frontSurfaceTensionForceModel

SourceFiles
    frontSurfaceTensionForceModel.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Interface for the front surface tension models. 

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


#include "frontSurfaceTensionForceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontSurfaceTensionForceModel, 0);
    defineRunTimeSelectionTable(frontSurfaceTensionForceModel, Dictionary)

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
//
frontSurfaceTensionForceModel::frontSurfaceTensionForceModel(const dictionary& configDict)
    :
        filterFieldName_(configDict.get<word>("filterField"))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<frontSurfaceTensionForceModel>
frontSurfaceTensionForceModel::New(const dictionary& configDict)
{
    const word name = configDict.get<word>("type");

    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = DictionaryConstructorTable(name);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "frontSurfaceTensionForceModel",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // Construct the model and return the autoPtr to the object.
    return tmp<frontSurfaceTensionForceModel> (ctorPtr(configDict));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
