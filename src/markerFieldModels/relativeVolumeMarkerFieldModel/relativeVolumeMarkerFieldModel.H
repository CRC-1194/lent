/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::relativeVolumeMarkerFieldModel

Description

SourceFiles
    relativeVolumeMarkerFieldModel.C

\*---------------------------------------------------------------------------*/

#ifndef relativeVolumeMarkerFieldModel_H
#define relativeVolumeMarkerFieldModel_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "sharpMarkerFieldModel.H"
// Let's first see if the header below is actually required
//#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class relativeVolumeMarkerFieldModel Declaration
\*---------------------------------------------------------------------------*/

class relativeVolumeMarkerFieldModel
:
    public sharpMarkerFieldModel
{

public:

    TypeName ("relativeVolume");

    // Constructors
    using markerFieldModel::markerFieldModel;

    // Member Functions
    virtual void calcMarkerField(volScalarField& markerField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
