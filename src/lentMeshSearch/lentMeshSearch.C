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

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#include "lentMeshSearch.H"
#include "dictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include <set>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(lentMeshSearch, 0); 
    defineRunTimeSelectionTable(lentMeshSearch, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentMeshSearch::lentMeshSearch(const dictionary& configDict)
{}

lentMeshSearch::lentMeshSearch()
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<lentMeshSearch>
lentMeshSearch::New(const dictionary& configDict)
{
    const word name = configDict.lookup("type"); 

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "lentMeshSearch::New(const word& name)"
        )   << "Unknown lentMeshSearch type "
            << name << nl << nl
            << "Valid lentMeshSearchs are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<lentMeshSearch> (cstrIter()(configDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lentMeshSearch::~lentMeshSearch()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label lentMeshSearch::cellContainingPoint(
    const point& p, 
    const fvMesh& mesh, 
    const label seedCell 
) const
{
    if (pointIsInCell(p, seedCell, mesh))
    {
        return seedCell;  
    }

    const volVectorField& C = mesh.C(); 
    scalar minDistance = mag(C[seedCell] - p);
    label minDistanceCell = seedCell; 

#if DEBUG
    Info << "minDistanceCell start = " << minDistanceCell << endl;
#endif

    // For all neighbour cells of the seed cell. 
    const labelListList& cellCells = mesh.cellCells(); 
    labelList neighborCells = cellCells[seedCell]; 

    const cellList& cells = mesh.cells(); 

    // Extend the tetrahedral stencil.
    if (cells[seedCell].size() == 4)
    {
        neighborCells = pointCellStencil(seedCell, mesh);  
    }

#if DEBUG
    Info << "used stencil = " << neighborCells  << endl;
#endif 

    forAll (neighborCells, I) 
    {
        label neighborCell = neighborCells[I]; 

#if DEBUG
        Info << "minDistance = " << minDistance << endl;
#endif

        if (pointIsInCell(p, neighborCell, mesh)) 
        {
            return neighborCell; 
        }

#if DEBUG
        Info << "checking neighbor cell = " << neighborCell << endl;
#endif
        scalar distance = mag(C[neighborCell] - p); 

#if DEBUG
        Info << "distance = " << distance << endl;
#endif 

        // Set label of the cell with the minimal distance. 
        if (distance <= minDistance) 
        {
            minDistance = distance; 
            minDistanceCell = neighborCell; 
#if DEBUG
            Info << "minDistanceCell = " << minDistanceCell << endl;
#endif 
        }
    }

    if (pointIsInCell(p, minDistanceCell, mesh)) 
    {
        return minDistanceCell; 
    } else 
    {

#if DEBUG
        Info << "skipping to cell " << minDistanceCell << endl;
#endif

        return cellContainingPoint(p, mesh, minDistanceCell); 
    }

    return -1; 
}

bool lentMeshSearch::pointIsInCell(
    const point p, 
    const label cellLabel, 
    const fvMesh& mesh, 
    scalar tolerance
) const
{
    const cellList& cells = mesh.cells(); 
    const cell& cell = cells[cellLabel];
    const labelList& own = mesh.faceOwner(); 
    const surfaceVectorField& Cf = mesh.Cf(); 
    const surfaceVectorField& Sf = mesh.Sf(); 

    bool pointIsInside = true;

#if DEBUG
    Info << "point " << p << endl;
#endif 

    // For all face labels of the cell.
    forAll (cell, I) 
    {
        label faceLabel = cell[I];

        vector faceNormal = Sf[faceLabel];   

        // If the cell does not own the face.
        if (! (cellLabel == own[cell[I]])) 
        {
            faceNormal *= -1;  
        }

        // Compute the vector from the face center to the point p.
        vector fp = p - Cf[cell[I]];  

        if ((fp & faceNormal) > tolerance) 
        {

#if DEBUG 
            Info << "point outside face = " << (fp & faceNormal) << endl; 
#endif
            pointIsInside = false;
            break;
        }

#if DEBUG
        else
        {
            Info << "point inside face = " << (fp & faceNormal) << endl; 
        }
#endif 

    }

    return pointIsInside;
}


labelList lentMeshSearch::pointCellStencil(
    label cellLabel, 
    const fvMesh& mesh
) const
{
    const faceList& faces = mesh.faces(); 
    const cellList& cells = mesh.cells(); 
    const labelListList& pointCells = mesh.pointCells(); 

    labelList cellPoints = cells[cellLabel].labels(faces); 

    std::set<label> newNeighborCells;

    forAll (cellPoints, I)
    {
        const labelList& addedNeighborCells = pointCells[cellPoints[I]]; 
        forAll(addedNeighborCells, J)
        {
            newNeighborCells.insert(addedNeighborCells[J]);
        }
    }

    return labelList(newNeighborCells.begin(), newNeighborCells.end());  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
