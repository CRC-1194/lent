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
    Foam::harmonicMarkerFieldModel

Description
    Abstract base class for the markerField function calculation from a signed 
    distance field.

SourceFiles
    harmonicMarkerFieldModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef harmonicMarkerFieldModel_H
#define harmonicMarkerFieldModel_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "markerFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class harmonicMarkerFieldModel Declaration
\*---------------------------------------------------------------------------*/

class harmonicMarkerFieldModel
    :
        public markerFieldModel
{

public:

    TypeName ("harmonic"); 

    // Constructors

        harmonicMarkerFieldModel(const dictionary& configDict);

    //- Destructor
    virtual ~harmonicMarkerFieldModel();

    // Member Functions
    
        virtual void calcMarkerField(
            volScalarField& markerField, 
            const volScalarField& signedDistance, 
            const volScalarField& searchDistanceSqr
        ) const;  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
