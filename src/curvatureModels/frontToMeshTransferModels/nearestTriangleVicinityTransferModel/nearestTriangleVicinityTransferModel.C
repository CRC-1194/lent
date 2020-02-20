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

// TODO:
// Use the functionality of the distance calculator and the octree search
// to set up the member functions of this class (TT)

#include "nearestTriangleVicinityTransferModel.H"

#include "addToRunTimeSelectionTable.H"
#include "fvcAverage.H"

#include "triSurfaceFrontGeoMesh.H"

#include <algorithm>
#include <assert.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(nearestTriangleVicinityTransferModel, 0);
    addToRunTimeSelectionTable(frontToMeshTransferModel, nearestTriangleVicinityTransferModel, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void nearestTriangleVicinityTransferModel::computeTrianglesInCellNeighbourhoodMapping(
    const fvMesh& mesh,
    const lentCommunication& communication
) const
{
    trianglesInCellNeighbourhood_.clear();

    const auto& triasInCell = communication.interfaceCellToTriangles();

    // FIXME: under which conditions can a cell contain a triangle and have
    // a markerfield value of 0 or 1 assuming all information is up-to-date? (TT)
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

        // TODO: if this assertion is violated, it means that there is bulk cell
        // (probably in the narrow band) that falsely assumes a value of
        // 0 < alpha < 1.
        // Using the interface cells to triangles / vertices connectivity could be
        // the basis (at least for the sharp markerfield model) to catch the
        // cells with erroneous alpha values
        assert(trianglesInNeighbourhood.size() > 0
                && "Error: Neighbourhood of interface cell contains no triangles");
    }
}

nearestTriangleVicinityTransferModel::labelPointPair nearestTriangleVicinityTransferModel::closestFrontElement(
    const point& P,
    const std::vector<label>& triangleIDs,
    const triSurfaceFront& front
) const
{
    // For now do not consider intersection to find the closest element
    // to a point. Instead use the triangle barycentres (TT)
    labelPointPair closestElement{-1, point{0.0, 0.0, 0.0}};

    const auto& triangleCentres = front.Cf();
    scalar minDistanceSquared = GREAT;

    for (const auto& triaID : triangleIDs)
    {
        if (magSqr(triangleCentres[triaID] - P) < minDistanceSquared)
        {
            minDistanceSquared = magSqr(triangleCentres[triaID] - P);
            closestElement.first = triaID;
            closestElement.second = triangleCentres[triaID];
        }
    }

    assert(closestElement.first >= 0 
            && "Error: no closest triangle found.");

    return closestElement;
}

label nearestTriangleVicinityTransferModel::cellContainingClosestElement(
    const label& triangleID,
    const lentCommunication& communication) const
{
    return communication.triangleToCell()[triangleID];
}

std::vector<label> nearestTriangleVicinityTransferModel::computeTrianglesInVicinity(
    const point& P,
    const std::vector<label>& candidates,
    const triSurfaceFront& front
) const
{
    std::vector<label> triangleIDs{};

    assert(candidates.size() > 0
            && "Error: list of candidates is empty");

    const auto searchRadiusSquared = computeSearchRadiusSquared();
    const auto& triangleCentres = front.Cf();

    for (const auto& triaID : candidates)
    {
        if (magSqr(triangleCentres[triaID] - P) <= searchRadiusSquared)
        {
            triangleIDs.push_back(triaID);
        }
    }

    assert(triangleIDs.size() > 0
            && "Error: empty triangle vicinity.");

    return triangleIDs;
}

scalar nearestTriangleVicinityTransferModel::weightedCurvatureAverage(
    const point& P,
    const std::vector<label>& triangleIDs,
    const triSurfaceFrontVectorField& curvatureNormals,
    const triSurfaceFront& front
) const
{
    scalar averageCurvature{0.0};

    const auto& triNormals = front.faceNormals();
    const auto& triCentres = front.Cf();
    scalar weightSum = 0.0;

    for (const auto& triaID : triangleIDs)
    {
        auto magCurvature = mag(curvatureNormals[triaID]);
        auto curvatureSign = sign(triNormals[triaID]&curvatureNormals[triaID]);
        auto triWeight = weight(mag(triCentres[triaID] - P));
        averageCurvature += triWeight*curvatureSign*magCurvature;
        weightSum += triWeight;
    }

    return averageCurvature/weightSum;
}

scalar nearestTriangleVicinityTransferModel::computeSearchRadiusSquared() const
{
    // Keep the option open to compute neighbourhood specific
    // search radii (TT)
    return searchRadiusSquared_;
}

scalar nearestTriangleVicinityTransferModel::weight(const scalar&) const
{
    // Simple average (TT)
    // TODO: Reasonable weighting approach with regard to distance? (TT)
    //
    // Normalize the distance with the search radius? Then the weighting problem
    // can be investigated on a unit sphere (TT)
    return 1.0;
}

void nearestTriangleVicinityTransferModel::setSearchRadiusSquared(const fvMesh& mesh) const
{
    // Assumption: the finest resolution is found in the narrow band around the
    // interface and the resolution in the narrow band is uniform.
    // Thus, take minimum delta as basis for the search radius computation
    if (searchRadiusSquared_ <= 0.0)
    {
        auto faceDeltasTmp = mesh.delta();
        auto minDelta = min(mag(faceDeltasTmp)).value();

        // TODO: the optimal choice of the search radius is an open question for
        // now. There are essentially two requirements:
        //  1) the search ball should lie in the bounds of the cell neighbourhood
        //  2) at least a single triangle must be located in the search ball
        //      to ensure that a curvature can be transfered
        //  (TT)
        searchRadiusSquared_ = pow(minDelta*searchRadiusCoefficient_, 2.0); 
    }
}

bool nearestTriangleVicinityTransferModel::isInterfaceCell(const scalar& volFraction) const
{
    return (volFraction > 0.0 && volFraction < 1.0);
}

// Taken from Dual Kriging interpolation
void nearestTriangleVicinityTransferModel::removeDuplicates(std::vector<label>& listOfLabels) const
{
    std::sort(listOfLabels.begin(), listOfLabels.end());
    std::vector<label>::iterator newEnd;
    newEnd = std::unique(listOfLabels.begin(), listOfLabels.end());
    listOfLabels.resize(std::distance(listOfLabels.begin(), newEnd));
}

std::vector<label> nearestTriangleVicinityTransferModel::cellNeighbourhood(const label& cellLabel, const fvMesh& mesh) const
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

nearestTriangleVicinityTransferModel::labelPair nearestTriangleVicinityTransferModel::patchIDAndLocalIndex(const label& faceID, const fvMesh& mesh) const
{
    const auto& faceCurvature = *curvatureBuffer(mesh);

    if (faceID < faceCurvature.size())
    {
        return labelPair{-1,faceID};
    }
    else
    {
        const auto& boundary = faceCurvature.boundaryField();

        forAll(boundary, patchI)
        {
            const auto& patch = boundary[patchI].patch();

            if (patch.start() <= faceID && faceID < patch.start() + patch.size())
            {
                return labelPair{patchI, faceID - patch.start()};
            }
        }
    }

    // Should never be reached, will cause a runtime error
    return labelPair{-1,-1};
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nearestTriangleVicinityTransferModel::nearestTriangleVicinityTransferModel(const dictionary& configDict)
:
    frontToMeshTransferModel{configDict},
    markerFieldName_{configDict.get<word>("markerFieldName")},
    searchRadiusCoefficient_{configDict.get<scalar>("searchRadiusCoefficient")},
    searchRadiusSquared_{-1.0}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void nearestTriangleVicinityTransferModel::transferCurvature(
    const triSurfaceFrontVectorField& curvatureNormals,
    const triSurfaceFront& front,
    const fvMesh& mesh
) const
{
    setBufferToZero(mesh);

    // Distribute curvature from front to Eulerian mesh
    auto& faceCurvature = *curvatureBuffer(mesh);

    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    );

    computeTrianglesInCellNeighbourhoodMapping(mesh, communication);
    setSearchRadiusSquared(mesh);
    
    // Statistics/debugging
    /*
    scalar averageNTrias = 0.0;

    for (const auto& cellTriasMap : trianglesInCellNeighbourhood_)
    {
        averageNTrias += cellTriasMap.second.size();
    }

    if (averageNTrias > 0.0)
    {
        averageNTrias /= trianglesInCellNeighbourhood_.size();
    }

    Info << nl << "# interface cells = " << trianglesInCellNeighbourhood_.size()
         << " with an average of " << averageNTrias
         << " triangles per neighbourhood" << nl << endl;
    */
    // end debugging
    
    const auto& faceCentres = mesh.Cf();
    const auto& cells = mesh.cells();
    const auto& markerField = mesh.lookupObject<volScalarField>(markerFieldName_);
    const auto& cellNearestTriangle = communication.cellsTriangleNearest();
    const auto& triangleToCell = communication.triangleToCell();

    forAll(markerField, cellID)
    {
        if (!isInterfaceCell(markerField[cellID]))
        {
            continue;
        }

        auto closestInterfaceCellID = triangleToCell[cellNearestTriangle[cellID].index()];
        {
            const auto& facesOfCell = cells[cellID];

            for (const auto& faceID : facesOfCell)
            {
                auto patchIDElementID = patchIDAndLocalIndex(faceID, mesh);

                // The following distiction is a consquence of how OpenFOAM
                // stores and provides access to boundary fields
                //
                // Inner faces
                if (patchIDElementID.first < 0)
                {
                    // Avoid duplicate curvature transfer for faces belonging to
                    // two interface cells
                    if (faceCurvature[faceID] != 0.0)
                    {
                        continue;
                    }

                    auto closestElementInfo = closestFrontElement(faceCentres[faceID], trianglesInCellNeighbourhood_[closestInterfaceCellID], front);
                    auto containingCellID = cellContainingClosestElement(closestElementInfo.first, communication);
                    auto trianglesInVicinity = computeTrianglesInVicinity(closestElementInfo.second, trianglesInCellNeighbourhood_.at(containingCellID), front);
                    faceCurvature[faceID] = weightedCurvatureAverage(closestElementInfo.second, trianglesInVicinity, curvatureNormals, front);
                }
                else
                // Boundary faces
                {
                    // Avoid duplicate curvature transfer for faces belonging to
                    // two interface cells
                    if (faceCurvature.boundaryField()[patchIDElementID.first][patchIDElementID.second] != 0.0)
                    {
                        continue;
                    }

                    const auto& faceCentre = faceCurvature.boundaryField()[patchIDElementID.first].patch().Cf()[patchIDElementID.second];

                    auto closestElementInfo = closestFrontElement(faceCentre, trianglesInCellNeighbourhood_[closestInterfaceCellID], front);
                    auto containingCellID = cellContainingClosestElement(closestElementInfo.first, communication);
                    auto trianglesInVicinity = computeTrianglesInVicinity(closestElementInfo.second, trianglesInCellNeighbourhood_.at(containingCellID), front);
                    faceCurvature.boundaryFieldRef()[patchIDElementID.first][patchIDElementID.second] = weightedCurvatureAverage(closestElementInfo.second, trianglesInVicinity, curvatureNormals, front);
                }
            }
        }
    }
}

std::shared_ptr<volScalarField> nearestTriangleVicinityTransferModel::cellCurvature(const fvMesh& mesh) const
{
    auto cellCurvatureTmp = fvc::average(*curvatureBuffer(mesh));

    return std::shared_ptr<volScalarField>{cellCurvatureTmp.ptr()};
}

std::shared_ptr<surfaceScalarField> nearestTriangleVicinityTransferModel::faceCurvature(const fvMesh& mesh) const
{
    return curvatureBuffer(mesh);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
