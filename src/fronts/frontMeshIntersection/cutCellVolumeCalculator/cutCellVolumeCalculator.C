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

scalar cutCellVolumeCalculator::cutCellVolume(const label& cellID) const
{
    // Note: this function computes the volume enclosed by the front
    // and the cell vertices with a negative signed distance
    scalar volume = 0.0;

    // Big steps:
    // 1) Find cell vertices with sd < 0
    // 2) Find faces completely located on negative side of front
    // 3) Reconstruct topology using edge vectors
    //      --> test if intesection point is on edge
    const volScalarField& signedDistance = 
        mesh_.lookupObject<volScalarField>(cellDistFieldName_);
    const pointScalarField& pointSignedDistance =
        mesh_.lookupObject<pointScalarField>(pointDistFieldName_);
    
    if (cellToTria_.found(cellID))
    {
        const faceList& faces = 
        pointField vertices(0);
        List<triFace> triangles(0);
        // Shit just got serious
        // ...
        //
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
