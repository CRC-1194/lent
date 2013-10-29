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

#if FULLDEBUG
    mutable cellSet visualizationCellSet_; 
    mutable label iterationCount_;
    void appendLabelAndWriteCellSet(label cellLabel) const; 
#endif

public:

    TypeName ("lentMeshSearch"); 

    declareRunTimeSelectionTable (
        tmp,  
        lentMeshSearch,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors
    
#if FULLDEBUG
        lentMeshSearch(const Time& runTime); 
#else
        lentMeshSearch();

        lentMeshSearch(const dictionary& configDict);
#endif 

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
