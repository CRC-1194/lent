/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 AUTHOR,AFFILIATION
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

#include "triangleInCellTransferModel.H"

#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"

#include "lentCommunication.H"
#include "triSurfaceFrontGeoMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(triangleInCellTransferModel, 0);
    addToRunTimeSelectionTable(frontToMeshTransferModel, triangleInCellTransferModel, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void triangleInCellTransferModel::initializeCellCurvatureField(const fvMesh& mesh) const
{
    if (cellCurvatureTmp_.empty())
    {
        const Time& runTime = mesh.time();  

        cellCurvatureTmp_ = 
            tmp<volScalarField> 
            (
                new volScalarField
                (
                    IOobject(
                        "cell_curvature", 
                        runTime.timeName(), 
                        mesh,
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    ), 
                    mesh, 
                    dimensionedScalar(
                        "zero", 
                        dimless/dimLength, 
                        0.0
                    )
                )
            );
    }
}

void triangleInCellTransferModel::resetCurvature() const
{
    auto& cellCurvature = cellCurvatureTmp_.ref();
    cellCurvature = dimensionedScalar{
                        "zero", 
                        dimless/dimLength, 
                        0.0
                    };
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triangleInCellTransferModel::triangleInCellTransferModel(const dictionary& configDict)
:
    frontToMeshTransferModel{configDict},
    cellCurvatureTmp_{}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void triangleInCellTransferModel::transferCurvature(
    const triSurfaceFrontVectorField& curvatureNormals,
    const triSurfaceFront& front,
    const fvMesh& mesh
) const
{
    initializeCellCurvatureField(mesh);
    resetCurvature();

    // Distribute curvature from front to Eulerian mesh
    auto& cellCurvature = cellCurvatureTmp_.ref();

    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 

    //------------------------------------------------------------------------
    // Arithmetic mean of all triangles in cell
    const auto& trianglesInCell = communication.interfaceCellToTriangles();
    const auto& faceNormal = front.Sf();

    for (const auto& cellTrianglesMap : trianglesInCell)
    {
        const auto& cellLabel = cellTrianglesMap.first;
        const auto& triangleLabels = cellTrianglesMap.second;

        for (const auto& tl : triangleLabels)
        {
            cellCurvature[cellLabel] += mag(curvatureNormals[tl])*sign(curvatureNormals[tl]&faceNormal[tl]);
        }

        cellCurvature[cellLabel] /= triangleLabels.size();
    }

    //------------------------------------------------------------------------
    // Propagate to non-interface cells
    const auto& cellToTriangle = communication.cellsTriangleNearest();
    const auto& triangleToCell = communication.triangleToCell();

    forAll(cellToTriangle, I)
    {
        const auto& hitObject = cellToTriangle[I];
        
        if (hitObject.hit())
        {
            cellCurvature[I] = cellCurvature[triangleToCell[hitObject.index()]];
        }
    }
}

tmp<volScalarField> triangleInCellTransferModel::cellCurvature() const
{
    const auto& cellCurvatureField = cellCurvatureTmp_.ref();
    return tmp<volScalarField>{cellCurvatureField};
}

tmp<surfaceScalarField> triangleInCellTransferModel::faceCurvature() const
{
    const auto& cellCurvatureField = cellCurvatureTmp_.ref();
    return fvc::interpolate(cellCurvatureField); 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
