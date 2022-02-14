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
    Foam::centerNormalConsistency

SourceFiles
    centerNormalConsistency.C

Author
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Makes the normals of the iso-surface consistent by using the signed
    distance values stored in cell centers.

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

#include "centerNormalConsistency.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "lentCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(centerNormalConsistency, 0);
    addToRunTimeSelectionTable(normalConsistency, centerNormalConsistency, Dictionary);

// * * * * * * * * * * * * * Private  member functions * * * * * * * * * * * //
void centerNormalConsistency::runNormalConsistencyAlgorithm(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField&
) const
{
    List<labelledTri>& triangles = static_cast<List<labelledTri>& > (front);
    const vectorField& triangleNormals = front.faceNormals();
    const pointField& frontPoints = front.points();
    const fvMesh& mesh = signedDistance.mesh();
    const pointField& cellCenters = mesh.C();

    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 
    const auto& triangleToCell = communication.triangleToCell();  

    // For all faces
    forAll (triangles, triangleI)
    {
        point P = frontPoints[triangles[triangleI][0]];
        label cellI = triangleToCell[triangleI];
        const point& C = cellCenters[cellI];
        const vector& n = triangleNormals[triangleI];

        if (((C - P) & n) * signedDistance[cellI] <  0)
        {
            triangles[triangleI].flip();
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

centerNormalConsistency::centerNormalConsistency(const dictionary& configDict)
    :
        normalConsistency(configDict)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
