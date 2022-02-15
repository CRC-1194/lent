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
    Foam::cutCellVolumeCalculator

SourceFiles
    cutCellVolumeCalculator.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Compute the phase volume in cells which are intersected by the front

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

#include "simpleTriangulator.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
point simpleTriangulator::geometricCentre(const labelList& pointIDs) const
{
    point centre(0.0, 0.0, 0.0);

    forAll(pointIDs, I)
    {
        centre += vertices_[pointIDs[I]];
    }

    centre /= pointIDs.size();

    return centre;
}

// Ensure bounds when usind std::acos
scalar simpleTriangulator::unitLimiter(scalar a) const
{
    if (a > 1.0) return 1.0;
    else if (a < -1.0) return -1.0;
    else return a;
}

scalar simpleTriangulator::angle(const vector& refEdge, const vector& a,
                                 scalar signedDistanceA) const
{
    scalar angle = std::acos(unitLimiter(refEdge & a / (mag(refEdge) * mag(a))));

    if (signedDistanceA < 0.0)
    {
        angle = 2*Foam::constant::mathematical::pi - angle;
    }

    return angle;
}

void simpleTriangulator::linkedSort(scalarList& reference,
                                    labelList& dependent) const
{
    // Sorts in ascending order
    // FIXME: inefficient sorting approach (TT)
    for (label I = 0; I < reference.size()-1; I++)
    {
        for (label K = I+1; K < reference.size(); K++)
        {
            if (reference[K] < reference[I])
            {
                std::swap(reference[I], reference[K]);
                std::swap(dependent[I], dependent[K]);
            }
        }
    }
}

void simpleTriangulator::orderPoints(labelList& pointIDs, const point& refPoint,
                                     const vector& normal) const
{
    vector refEdge = vertices_[pointIDs[0]] - refPoint;
    scalarList angles(pointIDs.size(), 0.0);
    analyticalPlane refPlane(refPoint, (normal ^ refEdge));

    forAll(pointIDs, I)
    {
        point& testPoint = vertices_[pointIDs[I]];
        vector testEdge = testPoint - refPoint;
        angles[I] = angle(refEdge, testEdge, refPlane.signedDistance(testPoint));
    }

    linkedSort(angles, pointIDs);
}

void simpleTriangulator::triangulate(labelList& pointIDs,
                                     const vector& faceNormal)
{
    point geoCentre = geometricCentre(pointIDs);
    vertices_.append(geoCentre);
    orderPoints(pointIDs, geoCentre, faceNormal);

    for (label I = 0; I < pointIDs.size()-1; I++)
    {
        triangles_.append
        (
            triFace(vertices_.size()-1, pointIDs[I], pointIDs[I+1])
        );
    }

    // Add last triangle manually to avoid list access violation
    triangles_.append
    (
        triFace(vertices_.size()-1, pointIDs.last(), pointIDs.first())
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
simpleTriangulator::simpleTriangulator(pointField& vertices,
                                       List<triFace>& triangles)
:
    vertices_(vertices),
    triangles_(triangles)
{
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
simpleTriangulator::~simpleTriangulator()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void simpleTriangulator::setVertices(pointField& vertices)
{
    vertices_ = vertices;
}

void simpleTriangulator::setTriangleList(List<triFace>& triangles)
{
    triangles_ = triangles;
}

void simpleTriangulator::triangulateFace(labelList& facePointIDs,
                                         const vector& faceNormal)
{
    triangulate(facePointIDs, faceNormal);
}

void simpleTriangulator::triangulateFace(labelList& facePointIDs, 
                                         const tmp<analyticalSurface>& surfaceTmp)
{
    point geoCentre = geometricCentre(facePointIDs);
    vector normal = surfaceTmp->normalToPoint(geoCentre);

    triangulate(facePointIDs, normal);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
