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

\*---------------------------------------------------------------------------*/

#include "FrontAndMeshCommunication.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mesh, class Front, class Connectivity>
Foam::frontTracking::FrontAndMeshCommunication<Mesh, Front, Connectivity>::FrontAndMeshCommunication (
    const Mesh& mesh, 
    const Front& front
)
:
    mesh_(mesh), 
    front_(front)
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//template<class Mesh, class Front, class Connectivity>
//Foam::frontTracking::FrontAndMeshCommunication<Mesh, Front, Connectivity>::~FrontAndMeshCommunication()
//{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Mesh, class Front, class Connectivity>
void
Foam::frontTracking::FrontAndMeshCommunication<Mesh, Front, Connectivity>::calcCellsToElementsDist(volScalarField& cellsToElementsDist)
{
    // Finds the minimal distance field between each cell and the front //

    // Get the map type
    typedef typename Connectivity::MapType MapType;

    // Get cells-elements connectivity
    const MapType& elementsToCells = connectivity_.elementsToCells(); 

    //MapType::const_iterator cellIt;  

    // For all cells-elements
    //for (cellIt = cellsToElements.begin(); cellIt != cellsToElements.end(); 
        //++cellIt)
    //{
        // Get cell-elements

        // Initialize min distance

        // For all cell-elements

            // Compute distance between cell centre and front element

            // If distance is lower than the minimal value
                // Minimal value becomes the current distance

        // Set the minimal distance between the cell and the front 
   //}
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//template<class Mesh, class Front, class Connectivity>
//void Foam::FrontAndMeshCommunication<Mesh, Front>::operator=(
    //const FrontAndMeshCommunication<Mesh, Front>& rhs
//)
//{
    //// Check for assignment to self
    //if (this == &rhs)
    //{
        //FatalErrorIn
        //(
            //"Foam::FrontAndMeshCommunication<Mesh, Front>::operator="
            //"(const Foam::FrontAndMeshCommunication<Mesh, Front>&)"
        //)   << "Attempted assignment to self"
            //<< abort(FatalError);
    //}
//}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
