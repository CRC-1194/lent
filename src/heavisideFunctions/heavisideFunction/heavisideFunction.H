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
    Foam::heavisideFunction

Description
    Abstract base class for the heaviside function calculation from a signed 
    distance field.

SourceFiles
    heavisideFunction.C

Authors:
    Tomislav Maric 
    tomislav.maric@gmx.com

\*---------------------------------------------------------------------------*/

#ifndef heavisideFunction_H
#define heavisideFunction_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heavisideFunction Declaration
\*---------------------------------------------------------------------------*/

class heavisideFunction
{

public:

    TypeName ("heavisideFunction"); 

    declareRunTimeSelectionTable (
        autoPtr,  
        heavisideFunction,  
        Empty,  
        (
            const word& name
        ),
        () 
    );

    // Constructors

        //- Construct null
        heavisideFunction();

    // Selectors

        //- Select null constructed
        static autoPtr<heavisideFunction> New(const word& name);


    //- Destructor
    virtual ~heavisideFunction();


    // Member Functions
        
        virtual void calcHeaviside(
            volScalarField& heaviside, 
            const volScalarField& signedDistance, 
            const volScalarField& narrowBandWidthSqr 
        ) const = 0;  

        bool distanceWithinNarrowBand(scalar distance, scalar narrowBandWidth) const; 

    // Member Operators 
    void operator=(const heavisideFunction& rhs); 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
