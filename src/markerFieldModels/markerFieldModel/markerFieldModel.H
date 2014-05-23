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
    Foam::markerFieldModel

Description
    Abstract base class for the markerField function calculation from a signed 
    distance field.

SourceFiles
    markerFieldModel.C

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#ifndef markerFieldModel_H
#define markerFieldModel_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class markerFieldModel Declaration
\*---------------------------------------------------------------------------*/

class markerFieldModel
    :
        public refCount
{

public:

    TypeName ("markerFieldModel"); 

    declareRunTimeSelectionTable (
        tmp,  
        markerFieldModel,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors

        markerFieldModel(const dictionary& configDict);

    // Selectors

        static tmp<markerFieldModel> New(const dictionary& configDict);

    // Destructor
        virtual ~markerFieldModel();

    // Member Functions
        
        virtual void calcMarkerField(
            volScalarField& markerField, 
            const volScalarField& signedDistance, 
            const volScalarField& searchDistanceSqr
        ) const = 0;  

        bool distanceWithinNarrowBand(scalar distance, scalar narrowBandWidth) const; 

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
