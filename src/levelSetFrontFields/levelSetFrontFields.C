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

\*---------------------------------------------------------------------------*/

#include "levelSetFrontFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace levelSetFrontTracking
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const word levelSetFrontLabelField::typeName("levelSetFrontLabelField");

template<>
const word levelSetFrontScalarField::typeName("levelSetFrontScalarField");

template<>
const word levelSetFrontVectorField::typeName("levelSetFrontVectorField");

template<>
const word levelSetFrontSphericalTensorField::typeName
("levelSetFrontSphericalTensorField");

template<>
const word levelSetFrontSymmTensorField::typeName
("levelSetFrontSymmTensorField");

template<>
const word levelSetFrontTensorField::typeName("levelSetFrontTensorField");


template<>
const word levelSetFrontPointLabelField::typeName("levelSetFrontPointLabelField");

template<>
const word levelSetFrontPointScalarField::typeName("levelSetFrontPointScalarField");

template<>
const word levelSetFrontPointVectorField::typeName("levelSetFrontPointVectorField");

template<>
const word levelSetFrontPointSphericalTensorField::typeName
("levelSetFrontPointSphericalTensorField");

template<>
const word levelSetFrontPointSymmTensorField::typeName
("levelSetFrontPointSymmTensorField");

template<>
const word levelSetFrontPointTensorField::typeName("levelSetFrontPointTensorField");


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking
} // End namespace Foam

// ************************************************************************* //
