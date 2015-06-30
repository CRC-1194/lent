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
    diffuseInterfaceProperties.C

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


#include "triSurfaceMeshDistanceFieldCalculator.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(triSurfaceMeshDistanceFieldCalculator, 0);
    addToRunTimeSelectionTable(
       distanceFieldCalculator,
       triSurfaceMeshDistanceFieldCalculator,
       Dictionary
    );

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceMeshDistanceFieldCalculator::triSurfaceMeshDistanceFieldCalculator(
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

void triSurfaceMeshDistanceFieldCalculator::calcCellsToFrontDistance(
    volScalarField& signedDistance,
    const volScalarField& searchDistanceSqr,
    const triSurfaceFront& front
)
{
    // FIXME: Used for parallelization, needs work. TM
    if (front.size() > 0)
    {
        signedDistance = dimensionedScalar(
            "GREAT",
            dimLength,
            GREAT
        );
    } else
    {
        signedDistance = dimensionedScalar(
            "-GREAT",
            dimLength,
            -GREAT
        );

    }

    const fvMesh& mesh = signedDistance.mesh();

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

    // FIXME: combine with the KVS algorithm. TM. 
    frontMesh.findNearest(
        C,
        searchDistanceSqr,
        cellsElementNearest_
    );

    // Create a list of the volume types: based on the cell centre, the
    List<volumeType> volType;
    // Fill the list of the volume types.
    frontMesh.getVolumeType(C, volType);

    // For all volume types.
    forAll(volType, I)
    {
        // Get the volume type.
        volumeType vT = volType[I];

        const pointIndexHit& h = cellsElementNearest_[I];

        if (h.hit())
        {
            // If the volume is OUTSIDE.
            if (vT == volumeType::OUTSIDE)
            {
                // Set the positive distance.
                signedDistance[I] = Foam::mag(C[I] - h.hitPoint());
            }
            // If the volume is inside.
            else if (vT == volumeType::INSIDE)
            {
                // Set the negative distance.
                signedDistance[I] = -Foam::mag(C[I] - h.hitPoint());
            }
        }
    }

    // TODO: for coupled / processor boundaries

    // TODO: GREAT --> make it a controllable configuraton value
    narrowBandTmp_->ensureNarrowBand(signedDistance, GREAT);
}

void triSurfaceMeshDistanceFieldCalculator::calcPointsToFrontDistance(
    pointScalarField& pointSignedDistance,
    const pointScalarField& pointSearchDistanceSqr,
    const triSurfaceFront& front
)
{
    pointSignedDistance.resize(pointSearchDistanceSqr.size());

    // FIXME: Used for parallelization, needs work. TM
    if (front.size() > 0)
    {
        pointSignedDistance = dimensionedScalar("GREAT", dimLength, GREAT);
    } else
    {
        pointSignedDistance = dimensionedScalar("-GREAT", dimLength, -GREAT);
    }

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

    List<volumeType> volType;
    frontMesh.getVolumeType(points, volType);

    // For all volume types.
    forAll(volType, I)
    {
        // Get the volume type.
        volumeType vT = volType[I];

        const pointIndexHit& h = pointsElementNearest_[I];

        if (h.hit())
        {
            // If the volume is OUTSIDE.
            if (vT == volumeType::OUTSIDE)
            {
                // Set the positive distance.
                pointSignedDistance[I] = Foam::mag(points[I] - h.hitPoint());
            }
            // If the volume is inside.
            else if (vT == volumeType::INSIDE)
            {
                // Set the negative distance.
                pointSignedDistance[I] = -Foam::mag(points[I] - h.hitPoint());
            }
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
