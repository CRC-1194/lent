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
    Foam::cutCellVolumeCalculator

SourceFiles
    cutCellVolumeCalculator.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Compute the phase volume in cells which are intersected by the front
    
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
        map(cellID) = List<label>(1, mappedID);
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
    assert (triaToCell.size() > 0);

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

// TODO: this function should be refactored (TT)
scalar cutCellVolumeCalculator::cutCellVolume(const label& cellID) const
{
    // Note: this function computes the volume enclosed by the front
    // and the cell vertices with a negative signed distance
    scalar volume = 0.0;

    const pointScalarField& pointSignedDistance =
        mesh_.lookupObject<pointScalarField>(pointDistFieldName_);
    const faceList& faces = mesh_.faces();
    const pointField& meshVertices = mesh_.points();

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
                                     provideNormal(cellFace, meshVertices));
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
                                      provideNormal(cellFace, meshVertices));
        }
    }

    volume = polyhedraVolume(vertices, triangles);

    assert (volume >= 0.0);
    assert (volume < mesh_.V()[cellID]);

    return volume;
}

scalar cutCellVolumeCalculator::tetVolume(const pointField& points,
                                          const triFace& base,
                                          const point& top) const
{
    return mag((points[base[0]] - top) &
               ((points[base[1]] - top) ^ (points[base[2]] - top)))/6.0;
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
                pointToPoint(T[K]) = fragmentVertices.size() - 1;
            }
        }

        fragmentTriangles.append
        (
            triFace(pointToPoint[T[0]], pointToPoint[T[1]], pointToPoint[T[2]])
        );
    }

    // Test mapping: first point of the the first triangle of each
    // triangle set has to be identical
    assert (refVertices[frontTriangles[frontTriangleIDs[0]][0]] 
            == fragmentVertices[fragmentTriangles[0][0]]);
}

label cutCellVolumeCalculator::facePosition(const face& cellFace,
                                    const pointScalarField& distance) const
{
    // Checks a face for three possible states:
    // 1) located on the positive side of the front (return 1)
    // 2) intersecting the front (return 0)
    // 3) located on the negative side of the front (return -1)
    bool hasPositiveVertex = false;
    bool hasNegativeVertex = false;
    
    forAll(cellFace, K)
    {
        if (distance[cellFace[K]] < 0.0)
        {
            hasNegativeVertex = true;
        }
        else
        {
            hasPositiveVertex = true;
        }
    }

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

bool cutCellVolumeCalculator::vectorsParallel(const vector& a, const vector& b) const
{
    if (fabs((a & b) / (mag(a)*mag(b)) - 1) < SMALL) return true;
    else return false;
}

vector cutCellVolumeCalculator::provideNormal(const face& cellFace,
                                              const pointField& vertices) const
{
    return (vertices[cellFace[1]] - vertices[cellFace[0]])
             ^ (vertices[cellFace[2]] - vertices[cellFace[0]]);
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
    return mesh_.V()[cellIndex] - cellVolumeNegativePhase(cellIndex);
}

scalar cutCellVolumeCalculator::cellVolumeNegativePhase(const label& cellIndex) const
{
    const volScalarField& signedDistance = 
        mesh_.lookupObject<volScalarField>(cellDistFieldName_);

    if (cellToTria_.found(cellIndex))
    {
        return cutCellVolume(cellIndex);
    }
    else
    {
        if (signedDistance[cellIndex] < 0.0)
        {
            return mesh_.V()[cellIndex];
        }
        else
        {
            return 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
