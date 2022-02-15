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
    Foam::frontToMeshTransferModel

SourceFiles
    frontToMeshTransferModel.C

Description
    Interface for different front to mesh transfer approaches intended
    for curvature propagation. 

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

#ifndef frontToMeshTransferModel_H
#define frontToMeshTransferModel_H


#include "typeInfo.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "triSurfaceFront.H"
#include "triSurfaceFrontFieldsFwd.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class frontToMeshTransferModel Declaration
\*---------------------------------------------------------------------------*/

class frontToMeshTransferModel
:
    public refCount
{

public:

    // Static data members
    TypeName("abstract");

    declareRunTimeSelectionTable (
        tmp,
        frontToMeshTransferModel,
        Dictionary,
        (
            const dictionary& configDict
        ),
        (configDict)
    )


    // Constructors
    frontToMeshTransferModel(const dictionary& configDict);


    // Selectors
    static tmp<frontToMeshTransferModel> New(const dictionary& configDict);


    //- Destructor
    virtual ~frontToMeshTransferModel() = default;


    // Member Functions
    virtual void transferCurvature(
        const triSurfaceFrontVectorField& curvatureNormals,
        const triSurfaceFront& front,
        const fvMesh& mesh
    ) const = 0;

    virtual std::shared_ptr<volScalarField> cellCurvature(const fvMesh&) const = 0;  

    virtual std::shared_ptr<surfaceScalarField> faceCurvature(const fvMesh&) const = 0; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
