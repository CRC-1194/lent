/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::optimizedOctreeDistanceCalculator

SourceFiles
    optimizedOctreeDistanceCalculator.H

Authors:
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Computes signed distance fields between the volume mesh and the immersed
    surface mesh using octree searches implemented in the optimizedOctree class.

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/


#ifndef optimizedOctreeDistanceCalculator_H
#define optimizedOctreeDistanceCalculator_H

#include "typeInfo.H"
#include "refCount.H"
#include "triSurfaceMesh.H"

#include "distanceFieldCalculator.H"
#include "narrowBandPropagation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class optimizedOctreeDistanceCalculator Declaration
\*---------------------------------------------------------------------------*/

class optimizedOctreeDistanceCalculator
    :
        public distanceFieldCalculator
{
    //- Narrow band propagator.
    tmp<narrowBandPropagation> narrowBandTmp_;

public:

    TypeName("optimizedOctree");

    // Constructors
        explicit optimizedOctreeDistanceCalculator(const dictionary& config);

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
