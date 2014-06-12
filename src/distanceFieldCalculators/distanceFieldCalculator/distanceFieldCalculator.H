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
    Foam::distanceFieldCalculator

Description

    Interface class for the distance field calculator. 

    Provides the interface for calculating the signed distance field between
    a volume mesh (cells and points) and an immersed surface mesh (elements).

SourceFiles
    distanceFieldCalculator.H

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef distanceFieldCalculator_H
#define distanceFieldCalculator_H

#include "typeInfo.H"
#include "refCount.H"

#include "triSurfaceFront.H"
#include "fvMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class distanceFieldCalculator Declaration
\*---------------------------------------------------------------------------*/

class distanceFieldCalculator
    :
        public refCount
{

    scalar narrowBandWidth_;

public:

    TypeName("distanceFieldCalculator"); 

    declareRunTimeSelectionTable (
        tmp,  
        distanceFieldCalculator,  
        Dictionary,  
        (
            const dictionary& configDict
        ),
        (configDict)
    );
     
    // Constructors
    
        explicit distanceFieldCalculator(
            const dictionary& configDict
        );

    // Destructor
    
        virtual ~distanceFieldCalculator();

    // Selector
    
        static tmp<distanceFieldCalculator> New(
            const dictionary& configDict
        );

    // Member Functions

        virtual void calcCellSearchDistance(
            volScalarField& searchDistanceSqr
        ); 

        virtual void calcPointSearchDistance(
            pointScalarField& pointSearchDistanceSqr,
            const volScalarField& cellSearchDistanceSqr
        ); 

        virtual void calcPointSearchDistance(
            pointScalarField& pointSearchDistanceSqr
        ); 

        virtual void calcCellsToFrontDistance(
            volScalarField& signedDistance, 
            const volScalarField& searchDistanceSqr, 
            const triSurfaceFront& front
        ) = 0;

        virtual void calcPointsToFrontDistance(
            pointScalarField& pointSignedDistance, 
            const pointScalarField& pointSearchDistanceSqr, 
            const triSurfaceFront& front
        ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
