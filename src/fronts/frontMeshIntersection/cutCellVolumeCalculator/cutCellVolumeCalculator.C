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

#include "cutCellVolumeCalculator.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void cutCellVolumeCalculator::addToMap(const label& cellID, 
                                       const label& mappedID,
                                       Map<List<label>>& map)
{
    if (map.found(cellID))
    {
        map[cellID].append(mappedID);
    }
    else
    {
        map[cellID] = List<label>(1, mappedID);
    }
}

// This method inverts the triangle-->cell mapping
void cutCellVolumeCalculator::cellToTriangle()
{
    const lentCommunication& communication =
        mesh_.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front_,mesh_)
    ); 

    const auto& triaToCell = communication.triangleToCell();

    forAll(triaToCell, I)
    {
        addToMap(triaToCell[I], I, cellToTria_);
    }
}

void cutCellVolumeCalculator::cellToFace()
{
    // Construct mapping cell-->face only for intesected cells
    List<label> intersectedCells = cellToTria_.toc();
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    forAll(intersectedCells, I)
    {
        forAll(owner, K)
        {
            if (owner[K] == intersectedCells[I])
            {
                addToMap(intersectedCells[I], K, cellToFace_);
            }
        }

        forAll(neighbour, K)
        {
            if (neighbour[K] == intersectedCells[I])
            {
                addToMap(intersectedCells[I], K, cellToFace_);
            }
        }
    }
}

// TODO: refactor this abomination of a function...
scalar cutCellVolumeCalculator::cutCellVolume(const label& cellID) const
{
    // Note: this function computes the volume enclosed by the front
    // and the cell vertices with a negative signed distance
    scalar volume = 0.0;

    const volScalarField& signedDistance = 
        mesh_.lookupObject<volScalarField>(cellDistFieldName_);
    const pointScalarField& pointSignedDistance =
        mesh_.lookupObject<pointScalarField>(pointDistFieldName_);
    const faceList& faces = mesh_.faces();
    const pointField& meshVertices = mesh_.points();
    const surfaceVectorField& faceNormals = mesh_.Sf();
    
    if (cellToTria_.found(cellID))
    {
        pointField vertices(0);
        List<triFace> triangles(0);

        const labelListList& faceToEdges = mesh_.faceEdges();
        const edgeList& edges = mesh_.edges();

        // Composition of the polyhedra is performed in three steps:
        // 1) Add the triangles of the front located in the given cell
        // 2) Triangulate and add faces of the cell which are completely
        //    located on the "negative" side (signed distance of all vertices
        //    < 0)
        // 3) Triangulate and add the intersected faces
        
        // Part 1
        frontFragment(cellToTria_[cellID], vertices, triangles);
        // It is necessary to be able to distinguish between front vertices
        // and cell vertices for the correct set up of the intersected faces
        label numFrontVertices = vertices.size();

        // Part 2 + 3
        labelList cellFaces = cellToFace_[cellID];
        simpleTriangulator triangulator(vertices, triangles);

        // NOTE: duplication of points in the following is intended as it
        // removes the necessity to establish a mapping between the gloabl
        // point list and the temporary local one
        forAll(cellFaces, I)
        {
            face cellFace = faces[cellFaces[I]];
            
            label facePos = facePosition(cellFace, pointSignedDistance);

            if (facePos == -1)
            {
                // triangulate complete face
                labelList pointIDs(cellFace.size());

                forAll(cellFace, K)
                {
                    vertices.append(meshVertices[cellFace[K]]);
                    pointIDs[K] = vertices.size() - 1;
                }

                triangulator.triangulateFace(pointIDs,
                                             faceNormals[cellFaces[I]]);
            }
            else if (facePos == 0)
            {
                // Set up face fraction and triangulate it
                labelList pointIDs(0);
                labelList faceEdges = faceToEdges[cellFaces[I]];

                // Add front vertices to point set which are located on
                // the current face
                forAll(faceEdges, K)
                {
                    edge E = edges[faceEdges[K]];

                    // intersected Edge
                    if (pointSignedDistance[E[0]]*pointSignedDistance[E[1]] < 0)
                    {
                       pointIDs.append(intersectionID(meshVertices[E[0]], 
                                        meshVertices[E[1]], vertices,
                                        numFrontVertices));
                    }
                }

                // Add all face vertices with negative signed distance
                forAll(cellFace, K)
                {
                    if (pointSignedDistance[cellFace[K]] < 0)
                    {
                        vertices.append(meshVertices[cellFace[K]]);
                        pointIDs.append(vertices.size() - 1);
                    }
                }
                
                triangulator.triangulateFace(pointIDs,
                                             faceNormals[cellFaces[I]]);
            }
        }

        volume = polyhedraVolume(vertices, triangles);
    }
    else
    {
        if (signedDistance[cellID] < 0.0)
        {
            volume = mesh_.V()[cellID];
        }
        else
        {
            volume = 0.0;
        }
    }

    return volume;
}

scalar cutCellVolumeCalculator::tetVolume(const pointField& points,
                                          const triFace& base,
                                          const point& top) const
{
    return mag((points[base[0]] - top) &
               ((points[base[1]] - top) ^ (points[base[2]])))/6.0;
}

scalar cutCellVolumeCalculator::polyhedraVolume(const pointField& points,
                                        const List<triFace>& triangles) const
{
    scalar volume = 0.0;

    point centre = geometricCentre(points);

    forAll(triangles, I)
    {
        volume += tetVolume(points, triangles[I], centre);
    }

    return volume;
}

point cutCellVolumeCalculator::geometricCentre(const pointField& points) const
{
    point centre(0.0, 0.0, 0.0);

    forAll(points, I)
    {
        centre += points[I];
    }

    return centre/points.size();
}

void cutCellVolumeCalculator::frontFragment(const labelList& frontTriangleIDs,
                                        pointField& fragmentVertices,
                                        List<triFace>& fragmentTriangles) const
{
    const pointField& refVertices = front_.localPoints();
    const List<labelledTri>& frontTriangles = front_.localFaces();
    Map<label> pointToPoint(frontTriangleIDs.size());

    // Set up mapping of global point list to local point list
    // and copy the triangles with the new point IDs
    forAll(frontTriangleIDs, I)
    {
        labelledTri T = frontTriangles[frontTriangleIDs[I]];

        forAll(T, K)
        {
            if (!pointToPoint.found(T[K]))
            {
                fragmentVertices.append(refVertices[T[K]]);
                pointToPoint[T[K]] = fragmentVertices.size() - 1;
            }
        }

        fragmentTriangles.append
        (
            triFace(pointToPoint[T[0]], pointToPoint[T[1]], pointToPoint[T[2]])
        );
    }
}

label cutCellVolumeCalculator::facePosition(const face& cellFace,
                                    const pointScalarField& distance) const
{
    // Checks a face for three possible states:
    // 1) located on the positive side of the front (return 1)
    // 2) intersecting the front (return 0)
    // 3 located on the negative side of the front (return -1)
    bool hasPositiveVertex = false;
    bool hasNegativeVertex = false;
    
    forAll(cellFace, K)
    {
        // TODO: add special treatment for points coinciding with the front.
        // Currently this case is not covered and funny things may happen
        if (distance[cellFace[K]] < 0.0)
        {
            hasNegativeVertex = true;
        }
        else
        {
            hasPositiveVertex = true;
        }
    }

    assert (hasPositiveVertex || hasNegativeVertex);

    if (hasPositiveVertex && !hasNegativeVertex) return 1;
    else if (hasPositiveVertex && hasNegativeVertex) return 0;
    else return -1;
}

label cutCellVolumeCalculator::intersectionID(const point& a, const point& b,
                         const pointField& vertices, const label limit) const
{
    label intersectID = -1;
    vector ab = b - a;
    vector aToVertex;

    for (label I = 0; I < limit; I++)
    {
        aToVertex = vertices[I] - a;
        
        if (vectorsParallel(ab, aToVertex))
        {
            intersectID = I;
            break;
        }
    }

    assert (intersectID > -1);
    return intersectID;
}

bool cutCellVolumeCalculator::vectorsParallel(const vector& a, vector b) const
{
    if ((a & b) / (mag(a)*mag(b)) - 1 < SMALL) return true;
    else return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
cutCellVolumeCalculator::cutCellVolumeCalculator(const fvMesh& mesh,
                                                 const triSurfaceFront& front,
                                                 const word& cellDistFieldName,
                                                 const word& pointDistFieldName)
:
    mesh_(mesh),
    front_(front),
    cellDistFieldName_(cellDistFieldName),
    pointDistFieldName_(pointDistFieldName),
    cellToTria_(front.localFaces().size()/3),
    cellToFace_(front.localFaces().size()/3)
{
    cellToTriangle();
    cellToFace();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
cutCellVolumeCalculator::~cutCellVolumeCalculator()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar cutCellVolumeCalculator::cellVolumePositivePhase(const label& cellIndex) const
{
    return mesh_.V()[cellIndex] - cutCellVolume(cellIndex);
}

scalar cutCellVolumeCalculator::cellVolumeNegativePhase(const label& cellIndex) const
{
    return cutCellVolume(cellIndex);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
