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
    Foam::triangleInCellTransferModel

Description

SourceFiles
    triangleInCellTransferModelI.H
    triangleInCellTransferModel.C
    triangleInCellTransferModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef triangleInCellTransferModel_H
#define triangleInCellTransferModel_H

#include "frontToMeshTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class triangleInCellTransferModel Declaration
\*---------------------------------------------------------------------------*/

class triangleInCellTransferModel
:
    public frontToMeshTransferModel
{
    // Private data
    mutable tmp<volScalarField> cellCurvatureTmp_;

    // Private Member Functions
    void initializeCellCurvatureField(const fvMesh&) const; 
    void resetCurvature() const;


public:

    // Static data members
    TypeName ("triangleInCell");


    // Constructors
    triangleInCellTransferModel(const dictionary& configDict);


    //- Destructor
    virtual ~triangleInCellTransferModel() = default;


    // Member Functions
    virtual void transferCurvature(
        const triSurfaceFrontVectorField& curvatureNormals,
        const triSurfaceFront& front,
        const fvMesh& mesh
    ) const;

    virtual tmp<volScalarField> cellCurvature() const;  

    virtual tmp<surfaceScalarField> faceCurvature() const; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
