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

#include "explicitDiffusionTransferModel.H"

#include "triSurfaceFrontGeoMesh.H"

#include "addToRunTimeSelectionTable.H"
#include "fvcAverage.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(explicitDiffusionTransferModel, 0);
    addToRunTimeSelectionTable(frontToMeshTransferModel, explicitDiffusionTransferModel, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// Duplication from nearestTriangleVicinityTransferModel
void explicitDiffusionTransferModel::computeTrianglesInCellNeighbourhoodMapping(
    const fvMesh& mesh,
    const lentCommunication& communication
) const
{
    trianglesInCellNeighbourhood_.clear();

    const auto& triasInCell = communication.interfaceCellToTriangles();

    for (const auto& cellTriaMap : triasInCell)
    {
        trianglesInCellNeighbourhood_[cellTriaMap.first] = std::vector<label>{};
    }

    for (auto& cellTriaMap: trianglesInCellNeighbourhood_)
    {
        auto neighbourCellIDs = cellNeighbourhood(cellTriaMap.first, mesh);
        auto& trianglesInNeighbourhood = trianglesInCellNeighbourhood_[cellTriaMap.first];
        for (const auto& cellID : neighbourCellIDs)
        {
            if (triasInCell.find(cellID) != triasInCell.end())
            {
                const auto& containedTriangles = triasInCell.at(cellID);

                for (const auto& triaID : containedTriangles)
                {
                    trianglesInNeighbourhood.push_back(triaID);
                }
            }
        }
    }
}

void explicitDiffusionTransferModel::computeBoundaryAndInterfaceValues(
    const triSurfaceFrontVectorField& curvatureNormals,
    const triSurfaceFront& front,
    const fvMesh& mesh
) const
{
    interfaceCellCurvature_.clear();
    boundaryCellCurvature_.clear();

    const auto& triNormals = front.faceNormals();

    scalar accumulatedCurvature{0.0};

    // Add interface cells
    for (const auto& cellTriaMap : trianglesInCellNeighbourhood_)
    {
        accumulatedCurvature = 0.0;

        for (const auto& triaID : cellTriaMap.second)
        {
            accumulatedCurvature += sign(triNormals[triaID]&curvatureNormals[triaID])
                                    *mag(curvatureNormals[triaID]);
        }

        interfaceCellCurvature_[cellTriaMap.first] = accumulatedCurvature / cellTriaMap.second.size();
    }

    // Add boundary cells of narrow band
    // Use the propagation cell -> closest triangle -> containing cell -> curvature
    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    );

    const auto& cellToTriangle = communication.cellsTriangleNearest();
    const auto& triangleToCell = communication.triangleToCell();

    const auto& signedDistance = mesh.lookupObject<volScalarField>(signedDistanceName_);

    auto magDistGradTmp = mag(fvc::grad(signedDistance));
    const auto& magDistGrad = magDistGradTmp.ref();

    forAll(magDistGrad, I)
    {
        // TODO: reasonable threshold for gradient of signed distance for
        // boundary cells of the narrow band? (TT)
        if (magDistGrad[I] > 1000.0 && mag(signedDistance[I]) < GREAT)
        {
            auto interfaceCellID = triangleToCell[cellToTriangle[I].index()];

            boundaryCellCurvature_[I] = interfaceCellCurvature_.at(interfaceCellID);
        }
    }
}

void explicitDiffusionTransferModel::setInterfaceCurvature(const fvMesh& mesh) const
{
    auto& cellCurvatureField = *curvatureBuffer(mesh);

    for (const auto& cellIDCurvature : interfaceCellCurvature_)
    {
        cellCurvatureField[cellIDCurvature.first] = cellIDCurvature.second;
    }
}

void explicitDiffusionTransferModel::setBoundaryCurvature(const fvMesh& mesh) const
{
    auto& cellCurvatureField = *curvatureBuffer(mesh);

    for (const auto& cellIDCurvature : boundaryCellCurvature_)
    {
        cellCurvatureField[cellIDCurvature.first] = cellIDCurvature.second;
    }
}

// Taken from Dual Kriging interpolation
void explicitDiffusionTransferModel::removeDuplicates(std::vector<label>& listOfLabels) const
{
    std::sort(listOfLabels.begin(), listOfLabels.end());
    std::vector<label>::iterator newEnd;
    newEnd = std::unique(listOfLabels.begin(), listOfLabels.end());
    listOfLabels.resize(std::distance(listOfLabels.begin(), newEnd));
}

std::vector<label> explicitDiffusionTransferModel::cellNeighbourhood(const label& cellLabel, const fvMesh& mesh) const
{
    std::vector<label> neighbourCellLabels{};

    const auto& cellPoints = mesh.cellPoints()[cellLabel];

    for (const auto& pointLabel : cellPoints)
    {
        const auto& connectedCells = mesh.pointCells()[pointLabel];

        for (const auto& aCellLabel : connectedCells)
        {
            neighbourCellLabels.push_back(aCellLabel);
        }
    }

    removeDuplicates(neighbourCellLabels);

    return neighbourCellLabels;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
explicitDiffusionTransferModel::explicitDiffusionTransferModel(const dictionary& configDict)
:
    frontToMeshTransferModel{configDict},
    diffusionIterations_{readLabel(configDict.lookup("diffusionIterations"))},
    signedDistanceName_{configDict.lookup("signedDistanceField")}
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void explicitDiffusionTransferModel::transferCurvature(
    const triSurfaceFrontVectorField& curvatureNormals,
    const triSurfaceFront& front,
    const fvMesh& mesh
) const
{
    setBufferToZero(mesh);

    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    );

    // Setup cellInNeighbourhood for all interface cell (-> all cells containing triangles)
    computeTrianglesInCellNeighbourhoodMapping(mesh, communication);

    // Compute average cell-centred curvature for each interface cell
    //      -> buffer this in a map
    computeBoundaryAndInterfaceValues(curvatureNormals, front, mesh);

    // Distribute curvature in narrow band using lentCommunication
    auto& cellCurvatureField = *curvatureBuffer(mesh);

    const auto& cellToTriangle = communication.cellsTriangleNearest();
    const auto& triangleToCell = communication.triangleToCell();

    for (const auto& iCellCurvature : interfaceCellCurvature_)
    {
        cellCurvatureField[iCellCurvature.first] = iCellCurvature.second;
    }

    forAll(cellToTriangle, I)
    {
        const auto& hitObject = cellToTriangle[I];
        
        if (hitObject.hit())
        {
            cellCurvatureField[I] = cellCurvatureField[triangleToCell[hitObject.index()]];
        }
    }

    // Diffusion loop
    for (int I = 0; I != diffusionIterations_; ++I)
    {
        cellCurvatureField = fvc::average(fvc::interpolate(cellCurvatureField));

        setBoundaryCurvature(mesh);
    }

    /*
    const auto& cellToTriangle = communication.cellsTriangleNearest();
    const auto& triangleToCell = communication.triangleToCell();

    for (int I = 0; I != diffusionIterations_; ++I)
    {
        for (const auto& iCellCurvature : interfaceCellCurvature_)
        {
            cellCurvatureField[iCellCurvature.first] = iCellCurvature.second;
        }

        forAll(cellToTriangle, I)
        {
            const auto& hitObject = cellToTriangle[I];
            
            if (hitObject.hit())
            {
                cellCurvatureField[I] = cellCurvatureField[triangleToCell[hitObject.index()]];
            }
        }

        cellCurvatureField = fvc::average(fvc::interpolate(cellCurvatureField));
    }

    for (int I = 0; I != diffusionIterations_; ++I)
    {
        cellCurvatureField = fvc::average(fvc::interpolate(cellCurvatureField));
    }
    */
}

std::shared_ptr<volScalarField> explicitDiffusionTransferModel::cellCurvature(const fvMesh& mesh) const
{
    return curvatureBuffer(mesh);
}

std::shared_ptr<surfaceScalarField> explicitDiffusionTransferModel::faceCurvature(const fvMesh& mesh) const
{
    auto cellCurvatureTmp = fvc::interpolate(*curvatureBuffer(mesh));

    return std::shared_ptr<surfaceScalarField>{cellCurvatureTmp.ptr()};
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
