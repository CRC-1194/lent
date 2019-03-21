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

#include "compactLevelSetCurvatureModel.H"

#include "addToRunTimeSelectionTable.H"

#include "lentCommunication.H"
#include "triSurfaceFields.H"

#include <algorithm>
#include <iterator>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(compactLevelSetCurvatureModel, 0);
    addToRunTimeSelectionTable(curvatureModel, compactLevelSetCurvatureModel, Dictionary);
    
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void compactLevelSetCurvatureModel::computeCurvature(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const volScalarField& curvatureInputField = 
        mesh.lookupObject<volScalarField>(inputFieldName()); 

    auto cellCurvatureTmp = levelSetCurvature(curvatureInputField);

    auto& cellCurvatureField = *curvatureBuffer(mesh);
    // TODO: performance optimization: avoid copying of the curvature field (TT)
    cellCurvatureField = cellCurvatureTmp.ref();

    // Apply corrections
    if (distanceCorrection_ == "sphere")
    {
        applySphereCorrection(mesh, front);
    }
    else if (distanceCorrection_ == "interpolation")
    {
        applyInterpolationCorrection(mesh, front);
    }

    // Write to temporary front field to communicate curvature back to
    // Eulerian mesh
    triSurfaceScalarField curvatureBuffer
    {
        IOobject(
            "curvatureBuffer",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        front,
        dimensionedScalar{
            "zero",
            dimless/dimLength,
            0
        }
    };

    const auto& communication = mesh.lookupObject<lentCommunication>(
                                    lentCommunication::registeredName(front, mesh)
                                ); 
    const auto& triangleToCell = communication.triangleToCell();

    forAll(triangleToCell, index)
    {
        curvatureBuffer[index] = cellCurvatureField[triangleToCell[index]];
    }

    // Propagate compact curvature to narrow band cells
    cellCurvatureField *= 0.0;

    const auto& cellToTriangle = communication.cellsTriangleNearest();

    forAll(cellToTriangle, I)
    {
        const auto& hitObject = cellToTriangle[I];
        
        if (hitObject.hit())
        {
            cellCurvatureField[I] = curvatureBuffer[hitObject.index()];
        }
    }
}

void compactLevelSetCurvatureModel::applySphereCorrection(
    const fvMesh& mesh,
    const triSurfaceFront& front
) const
{
    const auto& communication = mesh.lookupObject<lentCommunication>(
                                    lentCommunication::registeredName(front, mesh)
                                ); 
    const auto& signedDistance = inputField(mesh);
    auto& curvature = *curvatureBuffer(mesh);

    // Only used as a list of all cells intersected by the interface
    const auto& interfaceCells = communication.interfaceCellToTriangles();

    for (const auto& p : interfaceCells)
    {
        const auto& temp = curvature[p.first];
        // TODO: What happens here if the curvature is approximately zero?
        // According to IEEE 754 this should yield a curvature of zero as result (TT)
        curvature[p.first] = 2.0/(2.0/temp + signedDistance[p.first]);
    }
}

void compactLevelSetCurvatureModel::applyInterpolationCorrection(
    const fvMesh& mesh,
    const triSurfaceFront& front
) const
{
    const auto& communication = mesh.lookupObject<lentCommunication>(
                                    lentCommunication::registeredName(front, mesh)
                                ); 
    const auto& signedDistance = inputField(mesh);
    auto& curvature = *curvatureBuffer(mesh);

    // Only used as a list of all cells intersected by the interface
    const auto& interfaceCells = communication.interfaceCellToTriangles();
    const auto& cellsTriangleNearest = communication.cellsTriangleNearest();
    const auto& C = mesh.C();

    for (const auto& p : interfaceCells)
    {
        const auto& cellLabel = p.first;

        // Idea: find neighbour cell for which the connection of the
        // two cell centres deviates the least from the interface normal
        //
        // Note: it is not sufficient to look only at the cells connected by a
        // face
        scalar maxCosine = 0.0;
        label interpolationCellLabel = -1;
        auto normal = cellsTriangleNearest[cellLabel].hitPoint() - C[cellLabel];
        normal /= (mag(normal) + SMALL);

        // No need to correct if cell centre is located very close to the
        // interface
        if (mag(signedDistance[cellLabel]) > SMALL)
        {
            auto extendedNeighbourCells = findNeighbourCells(mesh, cellLabel);

            for (const auto& neighbourCellLabel : extendedNeighbourCells)
            {
                if (signedDistance[cellLabel]*signedDistance[neighbourCellLabel] > 0.0)
                {
                    continue;
                }

                auto cc = C[cellLabel] - C[neighbourCellLabel];
                cc /= mag(cc);

                auto cosine = mag(cc & normal);

                if (cosine > maxCosine)
                {
                    maxCosine = cosine;
                    interpolationCellLabel = neighbourCellLabel;
                }
            }

            curvature[cellLabel] = 
                (
                 curvature[cellLabel]*mag(signedDistance[interpolationCellLabel])
                 +
                 curvature[interpolationCellLabel]*mag(signedDistance[cellLabel])
                )
                /
                (mag(signedDistance[cellLabel]) + mag(signedDistance[interpolationCellLabel]));
        }
    }
}

// TODO: this probably needs some serious optimization (TT)
std::vector<label> compactLevelSetCurvatureModel::findNeighbourCells(
    const fvMesh& mesh,
    const label& cellLabel
) const
{
    std::vector<label> neighbourCells{};

    const auto& vertexLabels = mesh.cellPoints()[cellLabel];
    const auto& vertexToCell = mesh.pointCells();

    for (const auto& vertexLabel : vertexLabels)
    {
        const auto& connectedCells = vertexToCell[vertexLabel];

        for (const auto& connectedCellID : connectedCells)
        {
            if (connectedCellID != cellLabel)
            {
                neighbourCells.push_back(connectedCellID);
            }
        }
    }

    // Remove duplicate entries
    std::sort(neighbourCells.begin(), neighbourCells.end());
    
    std::vector<label>::iterator newEnd = 
                std::unique(neighbourCells.begin(), neighbourCells.end());

    neighbourCells.resize(std::distance(neighbourCells.begin(), newEnd));

    return neighbourCells;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
compactLevelSetCurvatureModel::compactLevelSetCurvatureModel(const dictionary& configDict)
:
    divergenceBasedCurvatureModel{configDict},
    distanceCorrection_{configDict.lookup("distanceCorrection")}
{
    if (
            distanceCorrection_ != "sphere"
            && distanceCorrection_ != "interpolation"
            && distanceCorrection_ != "off"
       )
    {
        FatalErrorIn (
            "frontCompactDivGradCurvatureModel::cellCurvature(front, mesh)"
        )   << "Unknown distance correction "
            << distanceCorrection_
            << " for the keyword 'distanceCorrection'.\n"
            << "Valid options are 'sphere', 'interpolation' or 'off'."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
