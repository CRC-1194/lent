/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::analyticalRandomizedCircle

Description

SourceFiles
    analyticalRandomizedCircleI.H
    analyticalRandomizedCircle.C
    analyticalRandomizedCircleIO.C

\*---------------------------------------------------------------------------*/

#ifndef analyticalRandomizedCircle_H
#define analyticalRandomizedCircle_H

#include "noiseGenerator.H"

#include "analyticalCircle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class analyticalRandomizedCircle Declaration
\*---------------------------------------------------------------------------*/

class analyticalRandomizedCircle
:
    public analyticalCircle
{
    // Private data
    noiseGenerator noiseGen_;

    // Private Member Functions


public:

    TypeName ("randomizedCircle");

    // Constructors
    analyticalRandomizedCircle() = default;
    analyticalRandomizedCircle(const dictionary& configDict);
    
    //- Destructor
    virtual ~analyticalRandomizedCircle() = default;


    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
