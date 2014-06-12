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


#include "centerNormalConsistency.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(centerNormalConsistency, 0);
    addToRunTimeSelectionTable(normalConsistency, centerNormalConsistency, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

centerNormalConsistency::centerNormalConsistency(const dictionary& configDict)
    :
        normalConsistency(configDict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void centerNormalConsistency::makeFrontNormalsConsistent(
    triSurface& front,
    const labelList& triangleCells,
    const volScalarField& signedDistance
) const
{
    List<labelledTri>& triangles = static_cast<List<labelledTri>& > (front);
    const vectorField& triangleNormals = front.faceNormals();
    const pointField& frontPoints = front.points();
    const fvMesh& mesh = signedDistance.mesh();
    const pointField& cellCenters = mesh.C();

    // For all faces
    forAll (triangles, E)
    {
        point P = frontPoints[triangles[E][0]];
        label cellI = triangleCells[E];
        point C = cellCenters[cellI];
        vector n = triangleNormals[E];
        n /= mag(n);

        scalar normalCellDistance = (C - P) & n;

        //if (((normalCellDistance < 0 && signedDistance[cellI] > 0)  ||
            //(normalCellDistance > 0 && signedDistance[cellI] < 0)))
        if (normalCellDistance * signedDistance[cellI] <  0)
        {
            triangles[E].flip();
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
