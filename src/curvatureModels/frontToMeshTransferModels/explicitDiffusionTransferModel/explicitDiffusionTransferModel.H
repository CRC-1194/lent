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

Class
    Foam::explicitDiffusionTransferModel

Description

SourceFiles
    explicitDiffusionTransferModelI.H
    explicitDiffusionTransferModel.C
    explicitDiffusionTransferModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef explicitDiffusionTransferModel_H
#define explicitDiffusionTransferModel_H

#include "frontToMeshTransferModel.H"

#include "CurvatureBufferLogic.H"
#include "lentCommunication.H"

#include <map>
#include <utility>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class explicitDiffusionTransferModel Declaration
\*---------------------------------------------------------------------------*/

class explicitDiffusionTransferModel
:
    public frontToMeshTransferModel,
    protected CurvatureBufferLogic<volScalarField, fvMesh, scalar>
{

private:

    // Private data
    label diffusionIterations_;
    word signedDistanceName_;
    
    mutable std::map<label, std::vector<label>>trianglesInCellNeighbourhood_{};
    mutable std::map<label, scalar> interfaceCellCurvature_{};
    mutable std::map<label, scalar> boundaryCellCurvature_{};

    // Private Member Functions
    void computeTrianglesInCellNeighbourhoodMapping(
        const fvMesh& mesh,
        const lentCommunication& communication
    ) const;

    void computeBoundaryAndInterfaceValues(
        const triSurfaceFrontVectorField& curvatureNormals,
        const triSurfaceFront& front,
        const fvMesh& mesh
    ) const;

    void setInterfaceCurvature(const fvMesh&) const;
    void setBoundaryCurvature(const fvMesh&) const;

    // Taken from Dual Kriging interpolation
    void removeDuplicates(std::vector<label>& listOfLabels) const;
    std::vector<label> cellNeighbourhood(const label& cellLabel, const fvMesh& mesh) const;


public:

    TypeName ("explicitDiffusion");

    // Constructors
    explicitDiffusionTransferModel(const dictionary& configDict);


    //- Destructor
    virtual ~explicitDiffusionTransferModel(){};


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
