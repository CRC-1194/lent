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

Class
    Foam::exactCurvatureModel

SourceFiles
    exactCurvatureModelI.H
    exactCurvatureModel.C
    exactCurvatureModelIO.C

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

    virtual std::shared_ptr<volVectorField> cellInterfaceNormals(
        const fvMesh&,
        const triSurfaceFront&
    ) const;

    virtual scalar curvatureAtPoint(const point& P) const = 0;
    virtual vector normalAtPoint(const point& P) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
