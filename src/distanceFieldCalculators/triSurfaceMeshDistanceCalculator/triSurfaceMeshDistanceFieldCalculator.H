/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
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
    Foam::triSurfaceMeshDistanceFieldCalculator

SourceFiles
    triSurfaceMeshDistanceFieldCalculator.H

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Computes signed distance fields between the volume mesh and the immersed
    surface mesh using octree searches implemented in the triSurfaceMesh class.

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

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
        explicit triSurfaceMeshDistanceFieldCalculator(const dictionary& config);

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
