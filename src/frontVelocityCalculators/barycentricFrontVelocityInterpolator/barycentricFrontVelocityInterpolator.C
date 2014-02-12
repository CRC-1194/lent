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

#include "barycentricFrontVelocityInterpolator.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(barycentricFrontVelocityInterpolator, 0); 

    addToRunTimeSelectionTable(
        frontVelocityCalculator,
        barycentricFrontVelocityInterpolator,
        Dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

barycentricFrontVelocityInterpolator::barycentricFrontVelocityInterpolator(const dictionary& configDict)
:
    frontVelocityCalculator(configDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void barycentricFrontVelocityInterpolator::calcFrontVelocity(
    triSurfaceFrontVectorField& frontVelocity, 
    const volVectorField& U,
    labelList& elementCells
) const
{
    const triSurfaceFront& front = frontVelocity.mesh(); 

    frontVelocity.resize(front.nPoints()); 
    //frontVelocity = dimensionedVector("zero", dimLength/dimTime, vector(0,0,0));

    interpolationCellPoint<vector> barycentric(U); 

    const List<labelledTri>& elements = front.localFaces(); 
    const pointField& vertices = front.points(); 

    const fvMesh& mesh = U.mesh(); 

    //const lentMeshSearch& searchAlg = getSearchAlgorithm();  

    forAll (elementCells, elementI)
    {
        const triFace& element = elements[elementI]; 

        forAll (element, vertexI)
        {
            const point& vertex = vertices[element[vertexI]];  

            label foundCell = -1; 

            if (!pointIsInCell(vertex, elementCells[elementI], mesh))
            {
                //elementCells[elementI] = searchAlg.cellContainingPoint(
                    //vertex, 
                    //mesh,
                    //elementCells[elementI]
                //); 

                // Initial code that allows the outflow of a front from a domain, 
                // which is important for parallelisation and inflow/outflow BCs.
                // TODO: profiling, as this introduced additional branching into 
                //       the part of the calculation called *very often*.  
                //       I'm counting on the CPU branching manager here. 
                
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

            if (foundCell > 0)
            {
                frontVelocity[element[vertexI]] = barycentric.interpolate(
                    vertices[element[vertexI]],  
                    elementCells[elementI]
                );
            }
        }
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

barycentricFrontVelocityInterpolator::~barycentricFrontVelocityInterpolator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //