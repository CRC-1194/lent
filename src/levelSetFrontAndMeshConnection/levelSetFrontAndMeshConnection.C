/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Author
    Tomislav Maric
    tomislav.maric@gmx.com

\*---------------------------------------------------------------------------*/

#include "levelSetFrontAndMeshConnection.H"

namespace Foam {
namespace levelSetFrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void levelSetFrontAndMeshConnection::resetCellsDistanceSqr()
{
    cellSearchDistSqr_ = 0;
}

void levelSetFrontAndMeshConnection::calcCellsDistanceSqr() 
{
    // Get the owner cells.
    const labelList& owner = mesh_.owner(); 
    // Get the neighbour cells.
    const labelList& neighbour = mesh_.neighbour(); 
    // Get the cell centres
    const volVectorField& C = mesh_.C();

    // For all faces.
    forAll (owner, faceI)
    {
        // Get the ownerI.
        label ownerI = owner[faceI];
        // Get the neighbourI. 
        label neighbourI = neighbour[faceI];

        // Compute the double squared distance. 
        vector distance (
            C[neighbourI] - C[ownerI]
        );

        // Thickness of the impacted cells.
        scalar distanceSquared = distance & distance;

        // Add the distance squared to the owner. 
        cellSearchDistSqr_[ownerI] += distanceSquared;
        // Add the distance squared to the neighbour.
        cellSearchDistSqr_[neighbourI] += distanceSquared;
    }

    // Get the mesh cells. 
    const cellList& cells = mesh_.cells();  
    forAll (cells, cellI)
    {
        // Average the square distance. 
        cellSearchDistSqr_[cellI] /= cells[cellI].size(); 
    }
}

void levelSetFrontAndMeshConnection::calcCellsToElements()
{
    front_.findNearest(mesh_.C(), cellSearchDistSqr_, cellsElementNearest_);
}

void levelSetFrontAndMeshConnection::updateCellsToElements()
{
    // Using the exisitng cellsToElementsMap_ search for the element in the
    // surrounding cells using a dendritic search. 
}

void levelSetFrontAndMeshConnection::calcElementsToCells()
{
    // Invert cellsToElements.
}

void levelSetFrontAndMeshConnection::updateElementsToCells()
{
}

void levelSetFrontAndMeshConnection::calcCellsElements()
{
    // If the mesh is moving or is experiencing topological changes.
    if (mesh_.changing() || mesh_.moving())
    {
        // Reset the search distance squared. 
        resetCellsDistanceSqr(); 
        // Calculate the search distance. 
        calcCellsDistanceSqr(); 
        // Calculate cells->elements connectivity.
    }

    // Compute the nearest front element to each cell. 
    calcCellsToElements(); 
    // Calculate elements->cells connectivity. 
    calcElementsToCells(); 
}

void levelSetFrontAndMeshConnection::updateCellsElements()
{
    updateCellsToElements(); 
    updateElementsToCells(); 
}

void levelSetFrontAndMeshConnection::initConnectivity()
{
    // TODO: check if mesh is moving or changing
    // Calculate the search distances.  
    calcCellsDistanceSqr(); 
    // Calculates the cells->elements connectivity. 
    calcCellsToElements(); 
    // Calculate the elements->cells connectivity.
    calcElementsToCells(); 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

levelSetFrontAndMeshConnection::levelSetFrontAndMeshConnection (
    const fvMesh& mesh, 
    const levelSetFront& front
)
    :
        mesh_(mesh), 
        front_(front), 
        //cellsToElementsMap_(mesh.cells().size()),
        //elementsToCellsMap_(front.size()), 
        cellSearchDistSqr_(
           IOobject (
               "cellSearchDistSqr", 
               mesh.time().timeName(), 
               mesh, 
               IOobject::READ_IF_PRESENT, 
               IOobject::NO_WRITE
           ), 
           mesh, 
           dimensionedScalar ( 
               "zero", 
               dimLength, 
               0
           )
        ), 
        cellsElementNearest_(mesh.cells().size())
{
    // Initialize the connectivity. 
    initConnectivity(); 
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void levelSetFrontAndMeshConnection::calcConnectivity()
{
    // Calculate cells<->elements connectivity.
    calcCellsElements(); 
}

void levelSetFrontAndMeshConnection::updateConnectivity()
{
    // Update cells<->elements connectivity.
    updateCellsElements();
}

// ************************************************************************* //

} // End namespace frontTracking 

// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
