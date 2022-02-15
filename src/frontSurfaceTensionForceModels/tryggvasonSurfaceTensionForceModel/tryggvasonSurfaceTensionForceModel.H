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
    Foam::tryggvasonSurfaceTensionForceModel

SourceFiles
    tryggvasonSurfaceTensionForceModel.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description

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


#ifndef tryggvasonSurfaceTensionForceModel_H
#define tryggvasonSurfaceTensionForceModel_H

#include "frontSurfaceTensionForceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class tryggvasonSurfaceTensionForceModel Declaration
\*---------------------------------------------------------------------------*/

class tryggvasonSurfaceTensionForceModel
    :
        public frontSurfaceTensionForceModel
{

public:

    TypeName ("Tryggvason");

    // Constructors 
    tryggvasonSurfaceTensionForceModel() = default;
    tryggvasonSurfaceTensionForceModel(const dictionary&);  

    // Destructor
    virtual ~tryggvasonSurfaceTensionForceModel() {};

    // Member Functions
    virtual tmp<surfaceScalarField> faceSurfaceTensionForce(
        const fvMesh& mesh, 
        const triSurfaceFront& frontMesh 
    ) const override; 

    virtual tmp<volVectorField> cellSurfaceTensionForce(
        const fvMesh& mesh,  
        const triSurfaceFront& frontMesh 
    ) const override;  

    virtual tmp<fvMatrix<vector>> surfaceTensionImplicitPart(
        const volVectorField&,
        const volScalarField&,
        const triSurfaceFront&
    ) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
