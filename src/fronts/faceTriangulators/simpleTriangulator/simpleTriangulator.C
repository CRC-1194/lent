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
    Foam::simpleTriangulator

SourceFiles
    simpleTriangulator.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Triangulates a given point set by decomposition using the geometric
    centre of the point set.
    
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
                                         tmp<analyticalSurface> surfaceTmp)
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
