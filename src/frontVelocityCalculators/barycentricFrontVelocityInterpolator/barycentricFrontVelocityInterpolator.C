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
    Foam::barycentricFrontVelocityInterpolator

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Interpolate front velocity using barycentric interpolation.

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
