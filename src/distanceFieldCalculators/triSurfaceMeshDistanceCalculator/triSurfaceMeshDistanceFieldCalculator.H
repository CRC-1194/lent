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
    Foam::triSurfaceMeshDistanceFieldCalculator

Description
    Computes signed distance fields between the volume mesh and the immersed
    surface mesh using octree searches implemented in the triSurfaceMesh class.


SourceFiles
    triSurfaceMeshDistanceFieldCalculator.H

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMeshDistanceFieldCalculator_H
#define triSurfaceMeshDistanceFieldCalculator_H

#include "typeInfo.H"
#include "refCount.H"
#include "triSurfaceMesh.H"

#include "distanceFieldCalculator.H"
#include "narrowBandPropagation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class triSurfaceMeshDistanceFieldCalculator Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMeshDistanceFieldCalculator
    :
        public distanceFieldCalculator
{

    //- List of indices to the nearest elements to cells.
    DynamicList<pointIndexHit> cellsElementNearest_; 

    //- List of indices to the nearest elements to points.
    DynamicList<pointIndexHit> pointsElementNearest_; 

    //- Narrow band propagator.
    tmp<narrowBandPropagation> narrowBandTmp_; 

public:

    TypeName("triSurfaceMeshCalculator"); 
     
    // Constructors
    
        explicit triSurfaceMeshDistanceFieldCalculator(
            const dictionary& config
        );

    // Destructor
    
        virtual ~triSurfaceMeshDistanceFieldCalculator();

    // Member Functions
        void calcCellsToFrontDistance(
            volScalarField& signedDistance, 
            const volScalarField& searchDistanceSqr,
            const triSurfaceFront& front
        );

        virtual void calcPointsToFrontDistance(
            pointScalarField& pointSignedDistance, 
            const pointScalarField& pointSearchDistanceSqr, 
            const triSurfaceFront& front
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
