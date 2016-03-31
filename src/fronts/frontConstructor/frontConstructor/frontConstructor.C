/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "frontConstructor.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
bool frontConstructor::differentSign(scalar a, scalar b) const
{
    if (a*b < 0.0) return true;
    else return false;
}

bool frontConstructor::notInList(label ID, labelList& list) const
{
    forAll(list, I)
    {
        if (list[I] == ID)
        {
            return false;
        }
    }

    return true;
}

label frontConstructor::pointIndex(label ID, const labelList& list) const
{
    label result = -1;

    forAll(list, I)
    {
        if (list[I] == ID) result = I;
    }

    return result;
}

void frontConstructor::cellIntersections(const labelList& cellEdges,
                                         labelList& intersections) const
{
    forAll(cellEdges, K)
    {
        // Remember: the position of an edge in the list
        // 'intersectedEdges_' corresponds to the position of the
        // intersecting point stored in 'frontVertices_'
        label pointID = pointIndex(cellEdges[K], intersectedEdges_);

        if (pointID > -1)
        {
            intersections.append(pointID);
        }
    }
}

void frontConstructor::createTriangles()
{
    findIntersectedEdges();
    computeIntersections();
    findIntersectedCells();

    const labelListList& cellToEdge = mesh_.cellEdges();
    labelList intersections(0);

    faceTriangulator triangulation_(frontVertices_, frontTriangles_,
                                    surfaceTmp_);

    forAll(intersectedCells_, I)
    {
        const labelList& cellEdges = cellToEdge[intersectedCells_[I]];
        intersections.clear();

        cellIntersections(cellEdges, intersections);
        triangulation_.triangulateFace(intersections);
        // This would be the place to capture the mapping between
        // triangles and the cell I
        // TODO: No mapping is created here as it would only be useful for
        // a first time step. Thus save line and simply use the
        // functionality of lentCommunication to etablish the mapping
    }
}

void frontConstructor::signedDistance()
{
    // Implement when you have figured out how to use this frickin' cursed
    // pointScalarField...
}

void frontConstructor::findIntersectedEdges()
{
    intersectedEdges_.clear();

    const edgeList& edges = mesh_.edges();
    const pointField& vertices = mesh_.points();

    forAll(edges, I)
    {
        edge E = edges[I];
        scalar distPointA = surfaceTmp_->signedDistance(vertices[E[0]]);
        scalar distPointB = surfaceTmp_->signedDistance(vertices[E[1]]);

        // Pure intersection
        if (differentSign(distPointA, distPointB))
        {
            intersectedEdges_.append(I);
        }
        // FIXME: this has to be solved in another way. Currently,
        // the intersection aka mesh point is assigned to multiple edges
        // resulting in a non unique mapping --> will result in degenerate
        // triangles
        else if (distPointA == 0.0 || distPointB == 0.0)
        {
            intersectedEdges_.append(I);
        }
    }
}

void frontConstructor::findIntersectedCells()
{
    intersectedCells_.clear();

    const labelListList& edgeToCell = mesh_.edgeCells();

    forAll(intersectedEdges_, I)
    {
        label edgeID = intersectedEdges_[I];

        labelList cells = edgeToCell[edgeID];

        forAll(cells, K)
        {
            // Avoid duplicate cell entries
            if (notInList(cells[K], intersectedCells_))
            {
                intersectedCells_.append(cells[K]);
            }
        }
    }
}

void frontConstructor::computeIntersections()
{
    frontVertices_.clear();
    frontVertices_.setSize(intersectedEdges_.size());

    const edgeList& edges = mesh_.edges();
    const pointField& vertices = mesh_.points();

    forAll(intersectedEdges_, I)
    {
        edge E = edges[intersectedEdges_[I]];

        frontVertices_[I] = surfaceTmp_->intersection(vertices[E[0]],
                                            vertices[E[1]]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontConstructor::frontConstructor(const tmp<analyticalSurface> surfaceTmp,
                                   const fvMesh& mesh)
:
    surfaceTmp_(surfaceTmp),
    mesh_(mesh)
{
    frontTriangles_ = List<triFace>(0);
    intersectedEdges_ = labelList(0);
    intersectedCells_ = labelList(0);
    frontVertices_ = pointField(0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

frontConstructor::~frontConstructor()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
