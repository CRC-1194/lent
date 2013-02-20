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

#ifndef levelSetFrontFieldsFwd_H
#define levelSetFrontFieldsFwd_H

#include "fieldTypes.H"
#include "levelSetFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace levelSetFrontTracking
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
class DimensionedField;

class levelSetFrontGeoMesh;

typedef Foam::DimensionedField<label, levelSetFrontGeoMesh>
    levelSetFrontLabelField;
typedef Foam::DimensionedField<scalar, levelSetFrontGeoMesh>
    levelSetFrontScalarField;
typedef Foam::DimensionedField<vector, levelSetFrontGeoMesh>
    levelSetFrontVectorField;
typedef Foam::DimensionedField<sphericalTensor, levelSetFrontGeoMesh>
    levelSetFrontSphericalTensorField;
typedef Foam::DimensionedField<symmTensor, levelSetFrontGeoMesh>
    levelSetFrontSymmTensorField;
typedef Foam::DimensionedField<tensor, levelSetFrontGeoMesh>
    levelSetFrontTensorField;

class levelSetFrontPointGeoMesh;

typedef Foam::DimensionedField<label, levelSetFrontPointGeoMesh>
    levelSetFrontPointLabelField;
typedef Foam::DimensionedField<scalar, levelSetFrontPointGeoMesh>
    levelSetFrontPointScalarField;
typedef Foam::DimensionedField<vector, levelSetFrontPointGeoMesh>
    levelSetFrontPointVectorField;
typedef Foam::DimensionedField<sphericalTensor, levelSetFrontPointGeoMesh>
    levelSetFrontPointSphericalTensorField;
typedef Foam::DimensionedField<symmTensor, levelSetFrontPointGeoMesh>
    levelSetFrontPointSymmTensorField;
typedef Foam::DimensionedField<tensor, levelSetFrontPointGeoMesh>
    levelSetFrontPointTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
