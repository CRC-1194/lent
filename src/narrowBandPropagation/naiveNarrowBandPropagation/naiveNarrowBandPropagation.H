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
    Foam::naiveNarrowBandPropagation

Description

SourceFiles
    naiveNarrowBandPropagation.C

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef naiveNarrowBandPropagation_H
#define naiveNarrowBandPropagation_H

#include "narrowBandPropagation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

/*---------------------------------------------------------------------------*\
                         Class naiveNarrowBandPropagation Declaration
\*---------------------------------------------------------------------------*/

class naiveNarrowBandPropagation
:
    public narrowBandPropagation

{

public:

    TypeName("naiveNarrowBand"); 

    // Constructors 
    
        naiveNarrowBandPropagation(const dictionary& configDict);

    // Destructor 
        
        virtual ~naiveNarrowBandPropagation(); 

    // Member Functions 

        void ensureNarrowBand(volScalarField& signedDistance, scalar L) const; 

        void ensureNarrowBand(pointScalarField& pointSignedDistance, scalar L) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
