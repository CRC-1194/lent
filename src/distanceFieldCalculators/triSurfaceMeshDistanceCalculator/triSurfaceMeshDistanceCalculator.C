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
    Foam::triSurfaceMeshDistanceCalculator

SourceFiles
    triSurfaceMeshDistanceCalculator.H

Authors:
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Computes signed distance fields between the volume mesh and the immersed
    surface mesh using octree searches implemented in the triSurfaceMesh class.

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


#include "triSurfaceMeshDistanceCalculator.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(triSurfaceMeshDistanceCalculator, 0);
    addToRunTimeSelectionTable(
       distanceFieldCalculator,
       triSurfaceMeshDistanceCalculator,
       Dictionary
    );

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceMeshDistanceCalculator::triSurfaceMeshDistanceCalculator(
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

void triSurfaceMeshDistanceCalculator::calcCellsToFrontDistance(
    volScalarField& signedDistance,
    const volScalarField& searchDistanceSqr,
    const triSurfaceFront& front
)
{
    // FIXME: Used for parallelization, needs work. TM
    //if (front.size() > 0)
    //{
        signedDistance = dimensionedScalar(
            "GREAT",
            dimLength,
            GREAT
        );
    //} else
    //{
        //signedDistance = dimensionedScalar(
            //"-GREAT",
            //dimLength,
            //-GREAT
        //);

    //}

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

void triSurfaceMeshDistanceCalculator::calcPointsToFrontDistance(
    pointScalarField& pointSignedDistance,
    const pointScalarField& pointSearchDistanceSqr,
    const triSurfaceFront& front
)
{
    pointSignedDistance.resize(pointSearchDistanceSqr.size());

    // FIXME: Used for parallelization, needs work. TM
    //if (front.size() > 0)
    //{
        pointSignedDistance = dimensionedScalar("GREAT", dimLength, GREAT);
    //} else
    //{
        //pointSignedDistance = dimensionedScalar("-GREAT", dimLength, -GREAT);
    //}

    // Get the cell centres.
    const pointMesh& pMesh = pointSignedDistance.mesh();

    const pointField& points = pMesh().points();

    // FIXME: work with octree directly. 
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
