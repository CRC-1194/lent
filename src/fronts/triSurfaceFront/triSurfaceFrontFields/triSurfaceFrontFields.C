/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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


Author
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "triSurfaceFrontFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const word triSurfaceFrontLabelField::typeName("triSurfaceFrontLabelField");

template<>
const word triSurfaceFrontScalarField::typeName("triSurfaceFrontScalarField");

template<>
const word triSurfaceFrontVectorField::typeName("triSurfaceFrontVectorField");

template<>
const word triSurfaceFrontSphericalTensorField::typeName
("triSurfaceFrontSphericalTensorField");

template<>
const word triSurfaceFrontSymmTensorField::typeName
("triSurfaceFrontSymmTensorField");

template<>
const word triSurfaceFrontTensorField::typeName("triSurfaceFrontTensorField");


template<>
const word triSurfaceFrontPointLabelField::typeName("triSurfaceFrontPointLabelField");

template<>
const word triSurfaceFrontPointScalarField::typeName("triSurfaceFrontPointScalarField");

template<>
const word triSurfaceFrontPointVectorField::typeName("triSurfaceFrontPointVectorField");

template<>
const word triSurfaceFrontPointSphericalTensorField::typeName
("triSurfaceFrontPointSphericalTensorField");

template<>
const word triSurfaceFrontPointSymmTensorField::typeName
("triSurfaceFrontPointSymmTensorField");

template<>
const word triSurfaceFrontPointTensorField::typeName("triSurfaceFrontPointTensorField");


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking
} // End namespace Foam

// ************************************************************************* //
