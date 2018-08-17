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
    Foam::divergenceBasedCurvatureModel

Description

SourceFiles
    divergenceBasedCurvatureModelI.H
    divergenceBasedCurvatureModel.C
    divergenceBasedCurvatureModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef divergenceBasedCurvatureModel_H
#define divergenceBasedCurvatureModel_H


#include "typeInfo.H"
#include "volFields.H"

#include "curvatureModel.H"
#include "CurvatureBufferLogic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class divergenceBasedCurvatureModel Declaration
\*---------------------------------------------------------------------------*/

class divergenceBasedCurvatureModel
:
    public curvatureModel,
    protected CurvatureBufferLogic<volScalarField, fvMesh, scalar>
{

private:

    // Private data
    word inputFieldName_;

    // Private Member Functions
    virtual void computeCurvature(const fvMesh&, const triSurfaceFront&) const;

protected:
    tmp<volScalarField> levelSetCurvature(const volScalarField& levelSetField) const;

public:

    TypeName ("divGrad")

    // Constructors
    divergenceBasedCurvatureModel(const dictionary& configDict);

    //- Destructor
    virtual ~divergenceBasedCurvatureModel() = default;


    // Member Functions
    virtual std::shared_ptr<volScalarField> cellCurvature(
        const fvMesh&, 
        const triSurfaceFront&
    ) const;  

    virtual std::shared_ptr<surfaceScalarField> faceCurvature(
        const fvMesh&, 
        const triSurfaceFront&
    ) const; 
    
    const word& inputFieldName() const;
    const volScalarField& inputField(const fvMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
