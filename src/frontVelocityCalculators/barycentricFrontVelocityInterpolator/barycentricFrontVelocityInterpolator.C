/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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
    Foam::barycentricFrontVelocityCalculator

SourceFiles
    frontVelocityCalculator.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Barycentric front velocity interpolation.

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"


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
    triSurfacePointVectorField& frontVelocity,
    const volVectorField& U,
    labelList& elementCells
) const
{
    const triSurface& front = frontVelocity.mesh();

    frontVelocity.resize(front.nPoints());

    frontVelocity = dimensionedVector("zero", dimVelocity, vector(0,0,0)); 

    interpolationCellPoint<vector> barycentric(U);

    // FIXME: Replace the barycentricVelocityInterpolator by a derived class of 
    // lentInterpolation. TM.

    const List<labelledTri>& elements = front.localFaces();
    const pointField& vertices = front.points();

    const fvMesh& mesh = U.mesh();

    // Why not triSurfacePointVectorField? 
    // FIXME: Remove the search. 
    // FIXME: Investigate theh triSurfacePointVectorField are the faces only moved, or the points? 
    forAll (elementCells, elementI) // FIXME: Remove the search, update element cells in the lentCommunication class. TM. 
    {
        const triFace& element = elements[elementI];

        forAll (element, vertexI)
        {
            label foundCell = -1;

            const point& vertex = vertices[element[vertexI]]; 

            if (!pointIsInCell(vertex, elementCells[elementI], mesh)) 
            {
                foundCell = cellContainingPoint(
                    vertex,
                    mesh,
                    elementCells[elementI]
                );

                if (foundCell > 0)
                {
                    elementCells[elementI] = foundCell;
                }
            }
            else
            {
                // FIXME: Investigate what happens if KVS doesn't locate the point. Move to 
                // lentCommunication anyway. TM.
                //FatalErrorIn("barycentricFrontVelocityInterpolator::calcFrontVelocity")
                    //<< "Element cell not found." << endl;
                foundCell = elementCells[elementI];
            }

            if (foundCell > 0)
            {
                // Barycentric interpolation happens here, separate the rest of the code!
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
