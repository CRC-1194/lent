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
#include "lentCommunication.H"

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

    // Get the lent communication structure from the registry. 
    // Get the non-const reference because the nearest maps are modified by 
    // the distance calculation. TM. 
    const fvMesh& mesh = signedDistance.mesh();
    lentCommunication& comm = const_cast<lentCommunication&>(
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
        )
    );

    // Get all data required for distance calculation.
    // TODO: Check if this is a bottleneck and remove the mesh constructor.  TM.
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
    const volVectorField& C = mesh.C();
    auto& cellsTriangleNearest = comm.cellsTriangleNearest(); 

    // TODO: Check this for bottlenecks. Is it possible to update the 
    // cellsTriangleNearest by the KVS algorithm? TM. 
    frontMesh.findNearest(
        C,
        searchDistanceSqr,
        cellsTriangleNearest 
    );

    const auto& faceNormals = front.faceNormals(); 
    
    // Calculate the distance to the cell center as the distance to the 
    // nearest front triangle. 
    forAll(cellsTriangleNearest, cellI)
    {
        const pointIndexHit& h = cellsTriangleNearest[cellI];

        if (h.hit())
        {
            const auto& normalVector = faceNormals[h.index()]; 
            auto distanceVector = C[cellI] - h.hitPoint(); 
            auto distSign = sign(normalVector & distanceVector); 
            signedDistance[cellI] = distSign * mag(distanceVector);  
        }
    }

    narrowBandTmp_->ensureNarrowBand(signedDistance, GREAT);
}

void optimizedOctreeDistanceCalculator::calcPointsToFrontDistance(
    pointScalarField& pointSignedDistance,
    const pointScalarField& pointSearchDistanceSqr,
    const triSurfaceFront& front
)
{
    pointSignedDistance = dimensionedScalar("GREAT", dimLength, GREAT);

    // Get the lent communication structure from the registry. 
    // Get the non-const reference because the nearest maps are modified by 
    // the distance calculation. TM. 
    const polyMesh& mesh = pointSignedDistance.mesh()();
    lentCommunication& comm = const_cast<lentCommunication&>(
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
        )
    );

    // Get the cell centres.
    const pointMesh& pMesh = pointSignedDistance.mesh();
    const pointField& points = pMesh().points();

    // Get all data required for distance calculation.
    // TODO: Check if this is a bottleneck and remove the 
    // mesh constructor.  TM.
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
    auto& pointsTriangleNearest = comm.pointsTriangleNearest(); 

    frontMesh.findNearest(
        points,
        pointSearchDistanceSqr,
        pointsTriangleNearest 
    );
    const auto& faceNormals = front.faceNormals(); 

    // Calculate the distance to the cell corner point as the distance 
    // to the nearest front triangle. 
    forAll(pointsTriangleNearest, pointI)
    {
        const pointIndexHit& h = pointsTriangleNearest[pointI];

        if (h.hit())
        {
            // Set the distance.
            const auto& normalVector = faceNormals[h.index()];
            auto distanceVector = points[pointI] - h.hitPoint(); 
            auto distSign = sign(normalVector & distanceVector); 
            pointSignedDistance[pointI] = distSign * mag(distanceVector);  
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
