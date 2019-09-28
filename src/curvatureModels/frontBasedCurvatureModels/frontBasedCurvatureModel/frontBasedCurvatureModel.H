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
    Foam::frontBasedCurvatureModel

Description

SourceFiles
    frontBasedCurvatureModelI.H
    frontBasedCurvatureModel.C
    frontBasedCurvatureModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef frontBasedCurvatureModel_H
#define frontBasedCurvatureModel_H

#include "CurvatureBufferLogic.H"
#include "curvatureModel.H"
#include "frontToMeshTransferModel.H"

#include "surfaceFields.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class frontBasedCurvatureModel Declaration
\*---------------------------------------------------------------------------*/

class frontBasedCurvatureModel
:
    public curvatureModel,
    protected CurvatureBufferLogic<triSurfaceFrontVectorField, triSurfaceFront, vector>
{

private:

    // Private data
    tmp<frontToMeshTransferModel> frontToMeshTmp_;


    // Private Member Functions
    void updateCurvatureBuffers(const fvMesh&, const triSurfaceFront&) const;
    virtual void computeCurvature(const fvMesh&, const triSurfaceFront&) const = 0;


public:

    TypeName ("frontCurvatureInterface");

    // Constructors
    frontBasedCurvatureModel(const dictionary& configDict);


    //- Destructor
    virtual ~frontBasedCurvatureModel(){};


    // Member Functions
    virtual std::shared_ptr<volScalarField> cellCurvature(
        const fvMesh&, 
        const triSurfaceFront&
    ) const;  

    virtual std::shared_ptr<surfaceScalarField> faceCurvature(
        const fvMesh&, 
        const triSurfaceFront&
    ) const; 

    virtual std::shared_ptr<triSurfaceFrontVectorField> frontCurvatureNormal(
        const fvMesh&,
        const triSurfaceFront&
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
