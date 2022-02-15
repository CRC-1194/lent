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
    Foam::nearestTriangleVicinityTransferModel.H

SourceFiles
    nearestTriangleVicinityTransferModel.C

Description

    Directly transfer the curvature from the front to the cell faces of interface
    cells by finding the closest triangle in a cell neighbourhood.
    The curvature is computed as a weighted average of all triangles contained in
    a ball centred around the barycentre of the closest triangle to the centre of a
    cell face

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)
 
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

#ifndef nearestTriangleVicinityTransferModel_H
#define nearestTriangleVicinityTransferModel_H

#include "frontToMeshTransferModel.H"

#include "CurvatureBufferLogic.H"
#include "lentCommunication.H"

#include <map>
#include <vector>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class nearestTriangleVicinityTransferModel Declaration
\*---------------------------------------------------------------------------*/

class nearestTriangleVicinityTransferModel
:
    public frontToMeshTransferModel,
    protected CurvatureBufferLogic<surfaceScalarField, fvMesh, scalar>
{
public:

    // Type aliases
    using labelPointPair = std::pair<label, point>;
    using labelPair = std::pair<label, label>;

private:

    // Private data
    const word markerFieldName_;
    scalar searchRadiusCoefficient_;
    mutable scalar searchRadiusSquared_;

    // Buffer for communication data
    mutable std::map<label, std::vector<label>> trianglesInCellNeighbourhood_;

    // Private Member Functions
    void computeTrianglesInCellNeighbourhoodMapping(
                        const fvMesh& mesh,
                        const lentCommunication& communication
                    ) const;

    labelPointPair closestFrontElement(
                        const point& P,
                        const std::vector<label>& triangleIDs,
                        const triSurfaceFront& front
                    ) const;

    label cellContainingClosestElement(
                        const label& triangleID,
                        const lentCommunication& communication
                    ) const;

    std::vector<label> computeTrianglesInVicinity(
                        const point& P,
                        const std::vector<label>& candidates,
                        const triSurfaceFront& front
                    ) const;

    scalar weightedCurvatureAverage(
                        const point& P,
                        const std::vector<label>& triangleIDs,
                        const triSurfaceFrontVectorField& curvatureNormals,
                        const triSurfaceFront& front
                    ) const;

    scalar computeSearchRadiusSquared() const;
    scalar weight(const scalar& distance) const;
    void setSearchRadiusSquared(const fvMesh& mesh) const;
    bool isInterfaceCell(const scalar& volFraction) const;

    // Taken from Dual Kriging interpolation
    void removeDuplicates(std::vector<label>& listOfLabels) const;
    std::vector<label> cellNeighbourhood(const label& cellLabel, const fvMesh& mesh) const;

    // FIXME: workaround for access to boundary faces. More efficient
    // approach possible? (TT)
    labelPair patchIDAndLocalIndex(const label& faceID, const fvMesh& mesh) const;


public:

    // Static data members
    TypeName ("nearestTriangleVicinity");


    // Constructors
    nearestTriangleVicinityTransferModel(const dictionary& configDict);


    //- Destructor
    virtual ~nearestTriangleVicinityTransferModel() = default;


    // Member Functions
    virtual void transferCurvature(
        const triSurfaceFrontVectorField& curvatureNormals,
        const triSurfaceFront& front,
        const fvMesh& mesh
    ) const;

    virtual std::shared_ptr<volScalarField> cellCurvature(const fvMesh&) const;  

    virtual std::shared_ptr<surfaceScalarField> faceCurvature(const fvMesh&) const; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
