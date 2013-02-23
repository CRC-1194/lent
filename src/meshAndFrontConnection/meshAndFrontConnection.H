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

Class
    Foam::meshAndFrontConnection

Description

Author
    Tomislav Maric
    tomislav.maric@gmx.com

SourceFiles
    meshAndFrontConnectionI.H
    meshAndFrontConnection.C
    meshAndFrontConnectionIO.C

\*---------------------------------------------------------------------------*/

#ifndef meshAndFrontConnection_H
#define meshAndFrontConnection_H

//#include "DynamicList.H"
//#include "DLList.H"
//#include "levelSetFront.H"
//#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace frontTracking {

/*---------------------------------------------------------------------------*\
                         Class meshAndFrontConnection Declaration
\*---------------------------------------------------------------------------*/

template<class Front>
class meshAndFrontConnection
{

private: 
    // Private data
    
    const fvMesh& mesh_; 
    const Front& front_;
    
    //- Cells<->Elements 
    //MapType cellsToElementsMap_; 
    //MapType elementsToCellsMap_; 

    //- Distance field for cells.
    //volScalarField cellSearchDistSqr_;

    //- List of indices to the nearest elements to cells.
    //DynamicList<pointIndexHit> cellsElementNearest_; 

    // Private member functions
    
    //- Connectivity calculators: cells<-->elements 
    //void resetCellsDistanceSqr(); 
    //void calcCellsDistanceSqr(); 
    //void calcCellsToElements(); 
    //void updateCellsToElements(); 
    //void calcElementsToCells();
    //void updateElementsToCells(); 

    //void calcCellsElements(); 
    //void updateCellsElements(); 

    //void initConnectivity(); 

public:

    // Constructors

        //- Construct from components
        meshAndFrontConnection(const fvMesh& mesh, const Front& front);

    // Member Functions
    
        //- Access
        const fvMesh& mesh()
        {
            return mesh_; 
        }

        const Front& front()
        {
            return front_; 
        }

        //const DynamicList<pointIndexHit>&  cellsToElementsNearest()
        //{
            //return cellsElementNearest_;  
        //}

        //bool validElement(const pointIndexHit& el)
        //{
            //return el.hit(); 
        //}

        //- Edit
        //- Calculate new connectivity 
        //void calcConnectivity();
        //- Update connectivity
        //void updateConnectivity(); 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "meshAndFrontConnection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
