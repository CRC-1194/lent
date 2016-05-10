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
    Foam::frontMeshSearch

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Abstract base class for supplemental mesh search operations.
    fvMesh supports:

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


#include "frontMeshSearch.H"
#include "dictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include <set>

//#include <sstream>
//#include <iomanip>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontMeshSearch, 0);
    defineRunTimeSelectionTable(frontMeshSearch, Dictionary);
    addToRunTimeSelectionTable(frontMeshSearch, frontMeshSearch, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//frontMeshSearch::frontMeshSearch(const Time& runTime)
//:
    //lastSeedCell_(-1),
    //visualizationCellSet_(
        //IOobject(
            //"frontMeshSearchCells",
            //runTime.timeName(),
            //runTime,
            //IOobject::NO_READ,
            //IOobject::AUTO_WRITE
        //)
    //),
    //iterationCount_(0),
//{}

frontMeshSearch::frontMeshSearch(const dictionary& configDict)
:
    lastDistance_(GREAT)
{}

frontMeshSearch::frontMeshSearch()
:
    lastDistance_(GREAT)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<frontMeshSearch>
frontMeshSearch::New(const dictionary& configDict)
{
    const word name = configDict.lookup("type");

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "frontMeshSearch::New(const word& name)"
        )   << "Unknown frontMeshSearch type "
            << name << nl << nl
            << "Valid frontMeshSearchs are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<frontMeshSearch> (cstrIter()(configDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

frontMeshSearch::~frontMeshSearch()
{}

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
//void frontMeshSearch::appendLabelAndWriteCellSet(label cellLabel) const
//{
    //++iterationCount_;

    //visualizationCellSet_.insert(cellLabel);

    //std::stringstream ss;

    //ss << "frontMeshSearchCellSet" << "-" << std::setw(20) << std::setfill('0')
        //<< iterationCount_;

    //visualizationCellSet_.rename(ss.str());
    //visualizationCellSet_.write();
//}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label frontMeshSearch::cellContainingPoint(
    const point& p,
    const fvMesh& mesh,
    const label seedCell
) const
{
    //appendLabelAndWriteCellSet(seedCell);

    if (pointIsInCell(p, seedCell, mesh, -SMALL)) // FIXME: Tolerance data member? TM.
    {
        //appendLabelAndWriteCellSet(seedCell);
        return seedCell;
    }

    const volVectorField& C = mesh.C();
    scalar minDistance = mag(C[seedCell] - p);
    label minDistanceCell = seedCell;

    //Info << "minDistanceCell start = " << minDistanceCell << endl;

    // For all neighbour cells of the seed cell.
    const labelListList& cellCells = mesh.cellCells();
    labelList neighborCells = cellCells[seedCell];

    const cellList& cells = mesh.cells();

    // Extend the tetrahedral stencil.
    if (cells[seedCell].size() == 4)
    {
        neighborCells = pointCellStencil(seedCell, mesh);
    }

    //Info << "used stencil = " << neighborCells  << endl;

    forAll (neighborCells, I)
    {
        label neighborCell = neighborCells[I];

        //Info << "minDistance = " << minDistance << endl;

        if (pointIsInCell(p, neighborCell, mesh, -SMALL)) // FIXME: Tolerance data member? TM.
        {
            //appendLabelAndWriteCellSet(neighborCell);
            lastDistance_ = minDistance;
            return neighborCell;
        }

        //Info << "checking neighbor cell = " << neighborCell << endl;
        scalar distance = mag(C[neighborCell] - p);

        //Info << "distance = " << distance << endl;

        // Set label of the cell with the minimal distance.
        if (distance <= minDistance)
        {
            minDistance = distance;
            minDistanceCell = neighborCell;
            //Info << "minDistanceCell = " << minDistanceCell << endl;
        }
    }

    if (pointIsInCell(p, minDistanceCell, mesh, -SMALL)) // FIXME: Tolerance data member? TM. 
    {
        //appendLabelAndWriteCellSet(seedCell);
        lastDistance_ = minDistance;
        return minDistanceCell;
    } else
    {
        if (mag(lastDistance_ - minDistance) < SMALL)
        {
            return minDistanceCell;
            //return -1;
        }
        else
        {
            //Info << "skipping to cell " << minDistanceCell << endl;
            lastDistance_ = minDistance;
            return cellContainingPoint(p, mesh, minDistanceCell); 
        }
    }

    return -1;
}

bool frontMeshSearch::pointIsInCell(
    const point p,
    const label cellLabel,
    const fvMesh& mesh,
    scalar tolerance
) const
{
    const cellList& cells = mesh.cells();
    const cell& cell = cells[cellLabel];
    const labelList& own = mesh.faceOwner();
    const vectorField& Cf = mesh.faceCentres();
    const vectorField& Sf = mesh.faceAreas();

    bool pointIsInside = true;

    //Info << "point " << p << endl;

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

            //Info << "point outside face = " << (fp & faceNormal) << endl;
            pointIsInside = false;
            break;
        }

        //else
        //{
            //Info << "point inside face = " << (fp & faceNormal) << endl;
        //}

    }

    return pointIsInside;
}

labelList frontMeshSearch::pointCellStencil(
    label cellLabel,
    const fvMesh& mesh
) const
{
    labelList result;

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
    // TODO: Improve efficiency. TM. 
    result.resize(newNeighborCells.size());

    std::copy(newNeighborCells.begin(), newNeighborCells.end(), result.begin()); 

    return result;
}

void frontMeshSearch::updateElementCells(
    DynamicList<label>& elementCells,
    const triSurfaceFront& front,
    const fvMesh& mesh
) const
{
    const List<labelledTri>& elements = front.localFaces();
    const pointField& vertices = front.points();

    forAll (elementCells, elementI)
    {
        const triFace& element = elements[elementI];

        forAll (element, vertexI)
        {
            label foundCell = -1;

            const point& vertex = vertices[element[vertexI]];

                if (!pointIsInCell(vertex, elementCells[elementI], mesh, -SMALL)) //FIXME: Tolerance data member? TM. 
            {
                foundCell  = cellContainingPoint(
                    vertex,
                    mesh,
                    elementCells[elementI]
                );

                if (foundCell > 0)
                {
                    elementCells[elementI] = foundCell;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
