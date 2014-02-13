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


Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#include "triSurfaceMeshDistanceFieldCalculator.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(triSurfaceMeshDistanceFieldCalculator, 0);
    addToRunTimeSelectionTable(
       lentDistanceFieldCalculator, 
       triSurfaceMeshDistanceFieldCalculator, 
       Dictionary
    );

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceMeshDistanceFieldCalculator::triSurfaceMeshDistanceFieldCalculator(
    const dictionary& config
)
:
    lentDistanceFieldCalculator(config),
    cellsElementNearest_(), 
    pointsElementNearest_(), 
    narrowBandTmp_( 
       narrowBandPropagation::New(config.subDict("narrowBandPropagation"))
    ) 
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triSurfaceMeshDistanceFieldCalculator::~triSurfaceMeshDistanceFieldCalculator()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurfaceMeshDistanceFieldCalculator::calcCellsToFrontDistance(
    volScalarField& signedDistance, 
    const volScalarField& searchDistanceSqr, 
    const triSurfaceFront& front
)
{
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

    frontMesh.findNearest(
        C, 
        searchDistanceSqr, 
        cellsElementNearest_
    );

    // Create a list of the volume types: based on the cell centre, the
    List<searchableSurface::volumeType> volType;
    // Fill the list of the volume types. 
    frontMesh.getVolumeType(C, volType);

    // For all volume types. 
    forAll(volType, I) 
    {
        // Get the volume type.
        searchableSurface::volumeType vT = volType[I];

        const pointIndexHit& h = cellsElementNearest_[I]; 

        if (h.hit()) 
        {
            // If the volume is OUTSIDE.
            if (vT == searchableSurface::OUTSIDE) 
            {
                // Set the positive distance.
                signedDistance[I] = Foam::mag(C[I] - h.hitPoint());
            }
            // If the volume is inside.
            else if (vT == searchableSurface::INSIDE) 
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

    frontMesh.findNearest(
        points, 
        pointSearchDistanceSqr, 
        pointsElementNearest_
    );

    List<searchableSurface::volumeType> volType;
    frontMesh.getVolumeType(points, volType);

    // For all volume types. 
    forAll(volType, I)
    {
        // Get the volume type.
        searchableSurface::volumeType vT = volType[I];

        const pointIndexHit& h = pointsElementNearest_[I]; 

        if (h.hit())
        {
            // If the volume is OUTSIDE.
            if (vT == searchableSurface::OUTSIDE)
            {
                // Set the positive distance.
                pointSignedDistance[I] = Foam::mag(points[I] - h.hitPoint());
            }
            // If the volume is inside.
            else if (vT == searchableSurface::INSIDE)
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
