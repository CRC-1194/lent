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
    Foam::optimizedOctreeDistanceCalculator

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Computes signed distance fields between the volume mesh and the immersed
    surface mesh using octree searches implemented in the optimizedOctree class.

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


#include "optimizedOctreeDistanceCalculator.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(optimizedOctreeDistanceCalculator, 0);
    addToRunTimeSelectionTable(
       distanceFieldCalculator,
       optimizedOctreeDistanceCalculator,
       Dictionary
    );

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

optimizedOctreeDistanceCalculator::optimizedOctreeDistanceCalculator(
    const dictionary& config
)
:
    distanceFieldCalculator(config),
    cellsElementNearest_(),
    pointsElementNearest_(),
    narrowBandTmp_(
       narrowBandPropagation::New(config.subDict("narrowBandPropagation"))
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void optimizedOctreeDistanceCalculator::calcCellsToFrontDistance(
    volScalarField& signedDistance,
    const volScalarField& searchDistanceSqr,
    const triSurfaceFront& front
)
{
    signedDistance = dimensionedScalar("GREAT", dimLength, GREAT);

    const fvMesh& mesh = signedDistance.mesh();

    // FIXME: Go away from the triSurfaceMesh copy, build and use the octree
    // based on the front directly. TM.
    triSurfaceMesh frontMesh(
        IOobject(
            "triSurfaceMesh",
            "frontMesh",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        front
    );

    // Get the cell centres.
    const volVectorField& C = mesh.C();

    // FIXME: Update the cellsElementNearest by the KVS algorithm. TM. 
    // Use the same list below to compute the distances from the point hit.
    frontMesh.findNearest(
        C,
        searchDistanceSqr,
        cellsElementNearest_
    );

    const auto& faceNormals = front.faceNormals(); 
    
    // For all cell-elements.  
    forAll(cellsElementNearest_, I)
    {
        const pointIndexHit& h = cellsElementNearest_[I];

        if (h.hit())
        {
            // Set the distance.
            const auto& normalVector = faceNormals[h.index()]; 
            auto distanceVector = C[I] - h.hitPoint(); 
            auto distSign = sign(normalVector & distanceVector); 
            signedDistance[I] = distSign * mag(distanceVector);  
        }
    }

    // TODO: GREAT --> make it a controllable configuraton value
    narrowBandTmp_->ensureNarrowBand(signedDistance, GREAT);
}

void optimizedOctreeDistanceCalculator::calcPointsToFrontDistance(
    pointScalarField& pointSignedDistance,
    const pointScalarField& pointSearchDistanceSqr,
    const triSurfaceFront& front
)
{
    pointSignedDistance.resize(pointSearchDistanceSqr.size());

    pointSignedDistance = dimensionedScalar("GREAT", dimLength, GREAT);

    // Get the cell centres.
    const pointMesh& pMesh = pointSignedDistance.mesh();

    const pointField& points = pMesh().points();

    triSurfaceMesh frontMesh(
        IOobject(
            "triSurfaceMesh",
            "frontMesh",
            pointSignedDistance.mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        front
    );

    // FIXME: combine with the KVS algorithm. TM. 
    frontMesh.findNearest(
        points,
        pointSearchDistanceSqr,
        pointsElementNearest_
    );

    // For all cell-elements.  
    const auto& faceNormals = front.faceNormals(); 
    forAll(pointsElementNearest_, I)
    {
        const pointIndexHit& h = pointsElementNearest_[I];

        if (h.hit())
        {
            // Set the distance.
            const auto& normalVector = faceNormals[h.index()];
            auto distanceVector = points[I] - h.hitPoint(); 
            auto distSign = sign(normalVector & distanceVector); 
            pointSignedDistance[I] = distSign * mag(distanceVector);  
        }
    }

    // TODO: GREAT --> make it a controllable configuraton value
    narrowBandTmp_->ensureNarrowBand(pointSignedDistance, GREAT);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
