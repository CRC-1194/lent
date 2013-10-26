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
    Foam::heavisideModel

Description
    Abstract base class for the heaviside function calculation from a signed 
    distance field.

SourceFiles
    heavisideModel.C

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#ifndef heavisideModel_H
#define heavisideModel_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class heavisideModel Declaration
\*---------------------------------------------------------------------------*/

class heavisideModel
    :
        public refCount
{

public:

    TypeName ("heavisideModel"); 

    declareRunTimeSelectionTable (
        tmp,  
        heavisideModel,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors

        heavisideModel(const dictionary& configDict);

    // Selectors

        static tmp<heavisideModel> New(const word& name, const dictionary& configDict);

    // Destructor
        virtual ~heavisideModel();

    // Member Functions
        
        virtual void calcHeaviside(
            volScalarField& heaviside, 
            const volScalarField& signedDistance
        ) const = 0;  

        bool distanceWithinNarrowBand(scalar distance, scalar narrowBandWidth) const; 

    // Member Operators 
    
        void operator=(const heavisideModel& rhs); 

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
