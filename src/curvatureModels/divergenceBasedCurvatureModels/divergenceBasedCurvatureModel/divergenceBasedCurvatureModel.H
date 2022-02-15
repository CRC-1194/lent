/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)
    Tomislav Maric (maric@mma.tu-darmstadt.de)
 
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
    ) const override;  

    virtual std::shared_ptr<surfaceScalarField> faceCurvature(
        const fvMesh&, 
        const triSurfaceFront&
    ) const override; 

    virtual std::shared_ptr<volVectorField> cellInterfaceNormals(
        const fvMesh&,
        const triSurfaceFront&
    ) const override;
    
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
