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

Authors
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

#include "triSurfaceFrontFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const word FrontTracking::triSurfaceFrontLabelField::typeName("triSurfaceFrontLabelField");

template<>
const word FrontTracking::triSurfaceFrontScalarField::typeName("triSurfaceFrontScalarField");

template<>
const word FrontTracking::triSurfaceFrontVectorField::typeName("triSurfaceFrontVectorField");

template<>
const word FrontTracking::triSurfaceFrontSphericalTensorField::typeName
("triSurfaceFrontSphericalTensorField");

template<>
const word FrontTracking::triSurfaceFrontSymmTensorField::typeName
("triSurfaceFrontSymmTensorField");

template<>
const word FrontTracking::triSurfaceFrontTensorField::typeName("triSurfaceFrontTensorField");


template<>
const word FrontTracking::triSurfaceFrontPointLabelField::typeName("triSurfaceFrontPointLabelField");

template<>
const word FrontTracking::triSurfaceFrontPointScalarField::typeName("triSurfaceFrontPointScalarField");

template<>
const word FrontTracking::triSurfaceFrontPointVectorField::typeName("triSurfaceFrontPointVectorField");

template<>
const word FrontTracking::triSurfaceFrontPointSphericalTensorField::typeName
("triSurfaceFrontPointSphericalTensorField");

template<>
const word FrontTracking::triSurfaceFrontPointSymmTensorField::typeName
("triSurfaceFrontPointSymmTensorField");

template<>
const word FrontTracking::triSurfaceFrontPointTensorField::typeName("triSurfaceFrontPointTensorField");


template<>
const word FrontTracking::triSurfaceFrontEdgeLabelField::typeName("triSurfaceFrontEdgeLabelField");

template<>
const word FrontTracking::triSurfaceFrontEdgeScalarField::typeName("triSurfaceFrontEdgeScalarField");

template<>
const word FrontTracking::triSurfaceFrontEdgeVectorField::typeName("triSurfaceFrontEdgeVectorField");

template<>
const word FrontTracking::triSurfaceFrontEdgeSphericalTensorField::typeName
("triSurfaceFrontEdgeSphericalTensorField");

template<>
const word FrontTracking::triSurfaceFrontEdgeSymmTensorField::typeName
("triSurfaceFrontEdgeSymmTensorField");

template<>
const word FrontTracking::triSurfaceFrontEdgeTensorField::typeName("triSurfaceFrontEdgeTensorField");


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
