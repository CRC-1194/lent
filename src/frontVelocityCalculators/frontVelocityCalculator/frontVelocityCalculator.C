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
    Foam::frontVelocityCalculator

Description
    Interface for the front velocity calculation. 

SourceFiles
    frontVelocityCalculator.C

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "frontVelocityCalculator.H"
#include "dictionary.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(frontVelocityCalculator, 0); 
    defineRunTimeSelectionTable(frontVelocityCalculator, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontVelocityCalculator::frontVelocityCalculator(const dictionary& configDict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<frontVelocityCalculator>
frontVelocityCalculator::New(const dictionary& configDict)
{
    const word name = configDict.lookup("type"); 

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "frontVelocityCalculator::New(const word& name)"
        )   << "Unknown frontVelocityCalculator type "
            << name << nl << nl
            << "Valid frontVelocityCalculators are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<frontVelocityCalculator> (cstrIter()(configDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

frontVelocityCalculator::~frontVelocityCalculator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//void frontVelocityCalculator::updateMeshCells(labelList& meshCells, const fvMesh& mesh)
//{
    //const triSurfaceFront& front = frontVelocity.mesh(); 

    //frontVelocity.resize(front.nPoints()); 

    //interpolationCellPoint<vector> barycentric(U); 

    //const List<labelledTri>& elements = front.localFaces(); 
    //const pointField& vertices = front.points(); 

    //const fvMesh& mesh = U.mesh(); 

    //forAll (elementCells, elementI)
    //{
        //const triFace& element = elements[elementI]; 

        //forAll (element, vertexI)
        //{
            //const point& vertex = vertices[element[vertexI]];  

            //label foundCell = -1; 

            //if (!pointIsInCell(vertex, elementCells[elementI], mesh))
            //{
                //foundCell  = cellContainingPoint(
                    //vertex, 
                    //mesh,
                    //elementCells[elementI]
                //); 

                //if (foundCell > 0)
                //{
                    //elementCells[elementI] = foundCell; 
                //}
            //}
        //}
    //}

//}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
