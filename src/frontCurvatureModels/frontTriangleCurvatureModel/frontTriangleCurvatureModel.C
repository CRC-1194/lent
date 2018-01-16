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

    if (curvatureNormalTmp_.empty())
    {
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
    else if (curvatureNormalTmp_->size() != front.UList<labelledTri>::size())
    {
        curvatureNormalTmp_->resize(front.UList<labelledTri>::size());
    }
}

void frontTriangleCurvatureModel::computeCurvature(const fvMesh& mesh, const triSurfaceFront& front) const
{
    initializeCurvatureNormal(mesh, front);

    auto& cn = curvatureNormalTmp_.ref();
    cn *= 0.0;

    const auto& normalCalculator = normalCalculatorTmp_.ref();
    auto frontVertexNormalsTmp = normalCalculator.vertexNormals(mesh, front);
    auto& n = frontVertexNormalsTmp.ref();

    const auto& faces = front.localFaces();
    const auto& p = front.localPoints();
    const auto& triArea = front.magSf();

    forAll(faces, I)
    {
        const auto& f = faces[I];

        cn[I] = 0.5*(
                        ((p[f[1]] - p[f[0]]) ^ (n[f[1]] + n[f[0]]))
                      + ((p[f[2]] - p[f[1]]) ^ (n[f[2]] + n[f[1]]))
                      + ((p[f[0]] - p[f[2]]) ^ (n[f[0]] + n[f[2]]))
                    )
                    / triArea[I];
    }

    auto& cellCurvature = cellCurvatureTmp_.ref();
    cellCurvature *= 0.0;

    // TODO: this kind of interpolation / transfer should be moved to
    // lentInterpolation or lentCommunication (TT)
    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 

    /*
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
    */

    // Test averaging
    const auto& trianglesInCell = communication.interfaceCellToTriangles();
    const auto& faceNormal = front.Sf();

    for (const auto& cellTrianglesMap : trianglesInCell)
    {
        const auto& cellLabel = cellTrianglesMap.first;
        const auto& triangleLabels = cellTrianglesMap.second;

        for (const auto& tl : triangleLabels)
        {
            cellCurvature[cellLabel] += mag(cn[tl])*sign(cn[tl]&faceNormal[tl]);
        }

        cellCurvature[cellLabel] /= triangleLabels.size();
    }

    // Propagate to non-interface cells
    const auto& cellToTriangle = communication.cellsTriangleNearest();
    const auto& triangleToCell = communication.triangleToCell();

    forAll(cellToTriangle, I)
    {
        const auto& hitObject = cellToTriangle[I];
        
        if (hitObject.hit())
        {
            cellCurvature[I] = cellCurvature[triangleToCell[hitObject.index()]];
            //cellCurvatureField[I] = curvatureBuffer[hitObject.index()];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
frontTriangleCurvatureModel::frontTriangleCurvatureModel(const dictionary& configDict)
:
    frontCurvatureModel{configDict},
    normalCalculatorTmp_{
        frontVertexNormalCalculator::New(configDict.subDict("normalCalculator"))
    },
    curvatureNormalTmp_{}
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
