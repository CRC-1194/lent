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
    Foam::sinusEdgeLengthNormalCalculator

SourceFiles
    sinusEdgeLengthNormalCalculator.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Compute the vertex normals as a weighted average of the normals of the
    connected triangles.
    The weights are computed as in interface tracking using as weight
    w = sin(alpha) / (|d0| * |d1|), where alpha is the angle at the vertex and
    d0 and d1 are the edges connecting the vertex with the two other points
    of the given triangle.

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

#include "sinusEdgeLengthNormalCalculator.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(sinusEdgeLengthNormalCalculator, 0);
    addToRunTimeSelectionTable(frontVertexNormalCalculator, sinusEdgeLengthNormalCalculator, Dictionary);

FixedList<label,2> sinusEdgeLengthNormalCalculator::notCentrePoint(const label& centreID, const face& trianglePointLabels) const
{
    FixedList<label,2> nonCentreLabels{};

    label I = 0;

    for (const auto& pointID : trianglePointLabels)
    {
        if (pointID != centreID)
        {
            nonCentreLabels[I] = pointID;
            ++I;
        }
    }

    return nonCentreLabels;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
sinusEdgeLengthNormalCalculator::sinusEdgeLengthNormalCalculator(const dictionary& configDict)
:
    frontVertexNormalCalculator{configDict}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<triSurfaceFrontPointVectorField> sinusEdgeLengthNormalCalculator::vertexNormals(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const Time& runTime = mesh.time();  

    tmp<triSurfaceFrontPointVectorField> normalsTmp
    (
        new triSurfaceFrontPointVectorField
        (
            IOobject(
                "frontNormals", 
                runTime.timeName(), 
                front,
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            front, 
            dimensionedVector(
                "zero", 
                dimless, 
                vector(0.0,0.0,0.0)
            )
        )
    );

    auto& normals = normalsTmp.ref();
    const auto& vertices = front.localPoints();
    const auto& faces = front.localFaces();
    const auto& pointToFaces = front.pointFaces();
    const auto& faceNormals = front.faceNormals();
    auto weight = 1.0;

    forAll(vertices, I)
    {
        const auto& v = vertices[I];
        const auto& connectedFaces = pointToFaces[I];

        for (const auto& fl : connectedFaces)
        {
            const auto& aFace = faces[fl];
            auto nonCentreLabels = notCentrePoint(I, aFace);

            auto d0 = vertices[nonCentreLabels[0]] - v;
            auto d1 = vertices[nonCentreLabels[1]] - v;
            
            weight = mag(d0^d1)/(magSqr(d0)*magSqr(d1));

            normals[I] += weight*faceNormals[fl];
        }
    }

    normals /= mag(normals);

    return normalsTmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
