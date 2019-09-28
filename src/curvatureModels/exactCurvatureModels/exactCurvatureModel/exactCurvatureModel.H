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
    Foam::exactCurvatureModel

Description

SourceFiles
    exactCurvatureModelI.H
    exactCurvatureModel.C
    exactCurvatureModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef exactCurvatureModel_H
#define exactCurvatureModel_H


#include "surfaceFields.H"
#include "typeInfo.H"

#include "CurvatureBufferLogic.H"
#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class exactCurvatureModel Declaration
\*---------------------------------------------------------------------------*/

class exactCurvatureModel
:
    public curvatureModel,
    protected CurvatureBufferLogic<surfaceScalarField, fvMesh, scalar>
{

private:

    // Private data
    const Switch write_; 
    mutable std::shared_ptr<volScalarField> cellCurvaturePtr_{nullptr};
    mutable bool cellCurvatureInitialized_{false};

    // Private Member Functions
    void initializeCellCurvatureBuffer(const fvMesh&) const;
    void updateCurvatureBuffers(const fvMesh&, const triSurfaceFront&) const;

public:

    TypeName ("exactCurvatureInterface")

    // Constructors
    exactCurvatureModel(const dictionary& configDict);

    
    // Destructor
    virtual ~exactCurvatureModel() = default;


    // Member Functions
    virtual std::shared_ptr<volScalarField> cellCurvature(
        const fvMesh&, 
        const triSurfaceFront&
    ) const;  

    virtual std::shared_ptr<surfaceScalarField> faceCurvature(
        const fvMesh&, 
        const triSurfaceFront&
    ) const; 

    virtual scalar curvatureAtPoint(const point& P) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
