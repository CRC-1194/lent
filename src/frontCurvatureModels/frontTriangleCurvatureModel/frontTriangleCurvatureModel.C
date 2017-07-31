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
    Foam::frontTriangleCurvatureModel

SourceFiles
    frontTriangleCurvatureModel.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Curvature model based on the surface tension model described in the
    2012 paper of Tukovic and Jasak

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

#include "frontTriangleCurvatureModel.H"

#include "addToRunTimeSelectionTable.H"

#include "lentCommunication.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontTriangleCurvatureModel, 0);
    addToRunTimeSelectionTable(frontCurvatureModel, frontTriangleCurvatureModel, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void frontTriangleCurvatureModel::initializeCurvatureNormal(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const Time& runTime = mesh.time();  

    curvatureNormalTmp_ = 
        tmp<triSurfaceFrontVectorField> 
        (
            new triSurfaceFrontVectorField
            (
                IOobject(
                    "curvature_normal", 
                    runTime.timeName(), 
                    front,
                    IOobject::NO_READ, 
                    IOobject::AUTO_WRITE
                ), 
                front, 
                dimensionedVector(
                    "zero", 
                    dimless, 
                    vector(0.0,0.0,0.0)
                )
            )
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
frontTriangleCurvatureModel::frontTriangleCurvatureModel(const dictionary& configDict)
:
    frontCurvatureModel{configDict},
    normalAlgorithm_{configDict.lookup("normals")},
    curvatureNormalTmp_{}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<triSurfaceFrontPointVectorField> frontTriangleCurvatureModel::sphereNormals(const fvMesh& mesh, const triSurfaceFront& front) const
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
    vector centre{0, 0, 0};

    const auto& points = front.localPoints();

    forAll(points, I)
    {
        centre += points[I];
    }
    centre /= points.size();

    forAll(points, I)
    {
        normals[I] = points[I] - centre;
    }
    normals /= mag(normals);

    return normalsTmp;
}

tmp<triSurfaceFrontPointVectorField> frontTriangleCurvatureModel::invDistanceNormals(const fvMesh& mesh, const triSurfaceFront& front) const
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
    const auto& centres = front.Cf();
    const auto& points = front.localPoints();
    const auto& pointToFaces = front.pointFaces();
    const auto& faceNormals = front.Sf();
    auto w = 1.0;

    forAll(points, I)
    {
        const auto& faces = pointToFaces[I];

        forAll(faces, K)
        {
            const auto& fl = faces[K];
            w = 1.0/mag(points[I] - centres[fl]);

            normals[I] += w*faceNormals[fl]/mag(faceNormals[fl]);
        }
    }

    normals /= mag(normals);

    return normalsTmp;
}

tmp<triSurfaceFrontPointVectorField> frontTriangleCurvatureModel::areaAveragedNormals(const fvMesh& mesh, const triSurfaceFront& front) const
{
    // For now: compute vertex normals as area weighted average of
    // normals of connected triangles
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

tmp<triSurfaceFrontPointVectorField> frontTriangleCurvatureModel::vertexNormals(const fvMesh& mesh, const triSurfaceFront&front) const
{
    if (normalAlgorithm_ == "areaAveraged")
    {
        return areaAveragedNormals(mesh, front);
    }
    else if (normalAlgorithm_ == "inversedDistance")
    {
        return invDistanceNormals(mesh, front);
    }
    else if (normalAlgorithm_ == "sphere")
    {
        return sphereNormals(mesh, front);
    }
    else
    {
        FatalErrorIn ("frontTriangleCurvatureModel::vertexNormals()")
            << "Unknown vertex normal algorithm '"
            << normalAlgorithm_ << "'\n"
            << "Valid options are: areaAveraged, inversedDistance, sphere"
            << exit(FatalError);
    }
}

tmp<volScalarField> frontTriangleCurvatureModel::cellCurvature(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const Time& runTime = mesh.time();  

    tmp<volScalarField> cellCurvatureTmp(
        new volScalarField(
            IOobject(
                "cellCurvature", 
                runTime.timeName(), 
                mesh, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh, 
            dimensionedScalar(
                "zero", 
                pow(dimLength, -1), 
                0
            )
        )
    );
    
    if (curvatureNormalTmp_.empty())
    {
        initializeCurvatureNormal(mesh, front);
    }

    auto& cn = curvatureNormalTmp_.ref();

    auto frontVertexNormalsTmp = vertexNormals(mesh, front);
    auto& n = frontVertexNormalsTmp.ref();

    const auto& faces = front.localFaces();
    const auto& p = front.localPoints();
    const auto& triArea = front.magSf();

    forAll(faces, I)
    {
        const auto& f = faces[I];
        cn[I] = 0.5*(
                        ((p[f[0]] - p[f[2]]) ^ (n[f[0]] + n[f[2]]))
                      + ((p[f[1]] - p[f[0]]) ^ (n[f[1]] + n[f[0]]))
                      + ((p[f[2]] - p[f[1]]) ^ (n[f[2]] + n[f[1]]))
                    )
                    / triArea[I];
    }

    auto& cellCurvature = cellCurvatureTmp.ref();

    // TODO: this kind of interpolation / transfer should be moved to
    // lentInterpolation or lentCommunication
    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 

    const auto& cellsTriangleNearest = communication.cellsTriangleNearest();
    const auto& faceNormal = front.Sf();

    forAll(cellsTriangleNearest, I)
    {
        const auto& hitObject = cellsTriangleNearest[I];

        if (hitObject.hit())
        {
            const label& fl = hitObject.index();

            cellCurvature[I] = mag(cn[fl])*sign(cn[fl]&faceNormal[fl]);
        }
    }

    return cellCurvatureTmp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
