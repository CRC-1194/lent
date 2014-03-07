/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::lentMeshSearch

Description
    Abstract base class for supplemental mesh search operations. 
    fvMesh supports:
    
        * linear search (for each cell, check if ...)
        * logarithmic search (uses octrees)

    For now, the lentMeshSearch algorithm is to implement variants of the 
    vicinity search algorithm that uses neighboring cell relationship to find a
    point in the mesh.

    Important Note: the commented out code is used for debugging purposes, as the
    class is extended with a cellset of found cells for debugging. The cell-set
    is updated within a recursive member function, so I see no way to extend the
    class for debugging using OOD principles. Let me know if a clean way of doing
    this exists that allows me to avoid modification to the class as done here.
    Tomislav
           
SourceFiles
    lentMeshSearch.C

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#ifndef lentMeshSearch_H
#define lentMeshSearch_H

#include "fvMesh.H"
#include "cellSet.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class lentMeshSearch Declaration
\*---------------------------------------------------------------------------*/

class lentMeshSearch
    :
        public refCount
{
    mutable scalar lastDistance_; 
    //mutable cellSet visualizationCellSet_; 
    //mutable label iterationCount_;
    //void appendLabelAndWriteCellSet(label cellLabel) const; 

public:

    TypeName ("lentMeshSearch"); 

    declareRunTimeSelectionTable(
        tmp,  
        lentMeshSearch,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors
    
        //lentMeshSearch(const Time& runTime); 
        
        lentMeshSearch();

        explicit lentMeshSearch(const dictionary& configDict);

    // Selectors

        static tmp<lentMeshSearch> New(const dictionary& configDict);

    // Destructor
        virtual ~lentMeshSearch();

    // Member Functions
        
        bool pointIsInCell(
            const point p, 
            const label cellLabel, 
            const fvMesh& mesh, 
            scalar tolerance=SMALL
        ) const; 

        virtual label cellContainingPoint(
            const point& p, 
            const fvMesh& mesh, 
            const label seedCell = -1
        ) const;  

        labelList pointCellStencil(
            label cellLabel, 
            const fvMesh& mesh
        ) const;

        // TM March 07 14
        // Updates elementCells using previous elementCells values as seed cells.  
       
        // Works directly for a static mesh - the lentMethodDynamic is made 
        // responsible for updating elementCells for topological changes of the
        // volume mesh.

        void updateElementCells(
            DynamicList<label>& elementCells, 
            const triSurfaceFront& front, 
            const fvMesh& mesh
        ) const; 

    // Member Operators 
    
        void operator=(const lentMeshSearch& rhs); 

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
