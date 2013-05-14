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

//#include "fvMeshAndFrontConnection.H"

namespace Foam {

namespace frontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


//void fvMeshAndFrontConnection::calcCellsElements()
//{
    //// If the mesh is moving or is experiencing topological changes.
    //if (mesh_.changing() || mesh_.moving())
    //{
        //// Reset the search distance squared. 
        //resetCellsDistanceSqr(); 
        //// Calculate the search distance. 
        //calcCellsDistanceSqr(); 
        //// Calculate cells->elements connectivity.
    //}

    //// Compute the nearest front element to each cell. 
    //calcCellsToElements(); 
    //// Calculate elements->cells connectivity. 
    //calcElementsToCells(); 
//}

//void fvMeshAndFrontConnection::updateCellsElements()
//{
    //updateCellsToElements(); 
    //updateElementsToCells(); 
//}

//void fvMeshAndFrontConnection::initConnectivity()
//{
    //// TODO: check if mesh is moving or changing
    //// Calculate the search distances.  
    //calcCellsDistanceSqr(); 
    //// Calculates the cells->elements connectivity. 
    //calcCellsToElements(); 
    //// Calculate the elements->cells connectivity.
    //calcElementsToCells(); 
//}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Front>
fvMeshAndFrontConnection<Front>::fvMeshAndFrontConnection (
    const fvMesh& mesh, 
    const Front& front
)
    :
        mesh_(mesh), 
        front_(front)
{}

template<class Front>
fvMeshAndFrontConnection<Front>::fvMeshAndFrontConnection (
    const Front& front,
    const fvMesh& mesh 
)
    :
        mesh_(mesh), 
        front_(front)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//void fvMeshAndFrontConnection::calcConnectivity()
//{
    //// Calculate cells<->elements connectivity.
    //calcCellsElements(); 
//}

//void fvMeshAndFrontConnection::updateConnectivity()
//{
    //// Update cells<->elements connectivity.
    //updateCellsElements();
//}

// ************************************************************************* //

} // End namespace frontTracking 

// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
