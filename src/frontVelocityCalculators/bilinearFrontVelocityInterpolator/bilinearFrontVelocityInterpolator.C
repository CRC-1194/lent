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

\*---------------------------------------------------------------------------*/

#include "bilinearFrontVelocityInterpolator.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(bilinearFrontVelocityInterpolator, 0); 

    addToRunTimeSelectionTable(
        frontVelocityCalculator,
        bilinearFrontVelocityInterpolator,
        Dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bilinearFrontVelocityInterpolator::bilinearFrontVelocityInterpolator(const dictionary& configDict)
:
    frontVelocityCalculator(configDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void bilinearFrontVelocityInterpolator::calcFrontVelocity(
    triSurfaceFrontVectorField& frontVelocity, 
    const volVectorField& U,
    labelList& elementCells
) const
{
    const triSurfaceFront& front = frontVelocity.mesh(); 

    frontVelocity.resize(front.nPoints()); 
    frontVelocity = dimensionedVector("zero", dimLength/dimTime, vector(0,0,0));

    interpolationCellPoint<vector> bilinear(U); 

    const List<labelledTri>& elements = front.localFaces(); 
    const pointField& vertices = front.points(); 

    const fvMesh& mesh = U.mesh(); 

    const lentMeshSearch& searchAlg = getSearchAlgorithm();  

    forAll (elementCells, elementI)
    {
        const triFace& element = elements[elementI]; 

        forAll (element, vertexI)
        {
            const point& vertex = vertices[element[vertexI]];  

            if (!searchAlg.pointInCellWithTolerance(vertex, elementCells[elementI], mesh))
            {
                elementCells[elementI] = searchAlg.cellContainingPoint(
                    vertex, 
                    mesh,
                    elementCells[elementI]
                ); 
            }

            frontVelocity[element[vertexI]] = bilinear.interpolate(
                vertices[element[vertexI]],  
                elementCells[elementI]
            );
        }
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

bilinearFrontVelocityInterpolator::~bilinearFrontVelocityInterpolator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
