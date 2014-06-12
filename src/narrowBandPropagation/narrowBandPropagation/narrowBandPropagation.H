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
    Foam::narrowBandPropagation

Description
    Class for propagating the switch in the sign of the sign distance field.

SourceFiles
    narrowBandPropagation.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef narrowBandPropagation_H
#define narrowBandPropagation_H

#include "refCount.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

/*---------------------------------------------------------------------------*\
                         Class narrowBandPropagation Declaration
\*---------------------------------------------------------------------------*/

class narrowBandPropagation
:
    public refCount
{

public:

    TypeName("narrowBandPropagation"); 

    declareRunTimeSelectionTable (
        tmp,  
        narrowBandPropagation,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors 
    
        narrowBandPropagation(const dictionary& configDict);

    // Destructor 
        virtual ~narrowBandPropagation(); 

    // Selectors

        static tmp<narrowBandPropagation> New(const dictionary& configDict);

    // Member Functions 

        virtual void ensureNarrowBand(
            volScalarField& signedDistance, 
            scalar L = GREAT
        ) const = 0; 

        virtual void ensureNarrowBand(
            pointScalarField& pointSignedDistance, 
            scalar L = GREAT
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
