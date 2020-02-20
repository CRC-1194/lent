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
    Foam::diffuseInterfaceProperties

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    Mathematical Modeling and Analysis
    Center of Smart Interfaces, TU Darmstadt

Description
    
    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#include "distanceNormalConsistency.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "lentCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(distanceNormalConsistency, 0);
    addToRunTimeSelectionTable(normalConsistency, distanceNormalConsistency, Dictionary);

// * * * * * * * * * * * * * Private  member functions * * * * * * * * * * * //
void distanceNormalConsistency::runNormalConsistencyAlgorithm(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField&
) const
{
    auto& triangles = static_cast<List<labelledTri>& > (front);
    const auto& triangleNormals = front.faceNormals();
    const auto& frontPoints = front.points();
    const auto& pointTriangles = front.pointFaces();
    const auto& mesh = signedDistance.mesh();
    const auto& meshPoints = mesh.points(); 
    const auto& cellCenters = mesh.C();

    const auto& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 

    // FIXME: Draft, refactor this and comment the code. TM. 

    const auto& cellsTriangleNearest = communication.cellsTriangleNearest(); 
    forAll (cellsTriangleNearest, cellI)
    {
        const auto triangleI = cellsTriangleNearest[cellI].index();

        if (triangleI > -1)
        {
            if (((cellCenters[cellI] - frontPoints[triangles[triangleI][0]]) &
                triangleNormals[triangleI]) < 0)
                triangles[triangleI].flip(); 

            forAll (triangles[triangleI], pointI)
            {
                const auto& pTriangles = pointTriangles[triangles[triangleI][pointI]]; 

                forAll (pTriangles, pTriangleI)
                {
                    if ((triangleNormals[pTriangles[pTriangleI]] 
                        & triangleNormals[triangleI]) < 0)
                        triangles[pTriangles[pTriangleI]].flip(); 
                }
            } 
        }
    }

    const auto& pointsTriangleNearest = communication.pointsTriangleNearest(); 
    forAll (pointsTriangleNearest, pointI)
    {
        const auto triangleI = pointsTriangleNearest[pointI].index();

        if (triangleI > -1)
        {
            if (((meshPoints[pointI] - frontPoints[triangles[triangleI][0]]) &
                triangleNormals[triangleI]) < 0)
                triangles[triangleI].flip(); 

            forAll (triangles[triangleI], pointI)
            {
                const auto& pTriangles = pointTriangles[triangles[triangleI][pointI]]; 

                forAll (pTriangles, pTriangleI)
                {
                    if ((triangleNormals[pTriangles[pTriangleI]] 
                        & triangleNormals[triangleI]) < 0)
                        triangles[pTriangles[pTriangleI]].flip(); 
                }
            } 
        }
    }


    //const auto& triangleToCell = communication.triangleToCell();  
    // For all faces
    //forAll (triangles, triangleI)
    //{
        //point P = frontPoints[triangles[triangleI][0]];
        //label cellI = triangleToCell[triangleI];
        //const point& C = cellCenters[cellI];
        //const vector& n = triangleNormals[triangleI];

        //if (((C - P) & n) * signedDistance[cellI] <  0)
        //{
            //triangles[triangleI].flip();
        //}
    //}
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

distanceNormalConsistency::distanceNormalConsistency(const dictionary& configDict)
    :
        normalConsistency(configDict)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
