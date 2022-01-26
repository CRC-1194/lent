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
    Foam::frontVertexNormalCalculator

SourceFiles
    frontVertexNormalCalculator.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Compute the normals at the front vertices as the area average of the
    surrounding triangles.

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

#include "frontVertexNormalCalculator.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontVertexNormalCalculator, 0);
    defineRunTimeSelectionTable(frontVertexNormalCalculator, Dictionary)
    addToRunTimeSelectionTable(frontVertexNormalCalculator, frontVertexNormalCalculator, Dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
frontVertexNormalCalculator::frontVertexNormalCalculator(const dictionary&)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
tmp<frontVertexNormalCalculator> frontVertexNormalCalculator::New(const dictionary& configDict)
{
    const word name = configDict.get<word>("type");

    auto* ctorPtr = DictionaryConstructorTable(name);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "frontVertexNormalCalculator",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return tmp<frontVertexNormalCalculator>(ctorPtr(configDict));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<triSurfaceFrontPointVectorField> frontVertexNormalCalculator::vertexNormals(const fvMesh& mesh, const triSurfaceFront& front) const
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
    const auto& faceNormals = front.Sf();
    const auto& faces = front.localFaces();

    forAll(faces, I)
    {
        const auto& aFace = faces[I];

        forAll(aFace, K)
        {
            normals[aFace[K]] +=faceNormals[I];
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
