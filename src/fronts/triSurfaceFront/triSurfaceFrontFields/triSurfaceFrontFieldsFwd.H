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


Authors
    Tomislav Maric maric@csi.tu-darmstadt.de
\*---------------------------------------------------------------------------*/

#ifndef triSurfaceFrontFieldsFwd_H
#define triSurfaceFrontFieldsFwd_H

#include "fieldTypes.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
class DimensionedField;

class triSurfaceFrontGeoMesh;

typedef Foam::DimensionedField<label, triSurfaceFrontGeoMesh>
    triSurfaceFrontLabelField;
typedef Foam::DimensionedField<scalar, triSurfaceFrontGeoMesh>
    triSurfaceFrontScalarField;
typedef Foam::DimensionedField<vector, triSurfaceFrontGeoMesh>
    triSurfaceFrontVectorField;
typedef Foam::DimensionedField<sphericalTensor, triSurfaceFrontGeoMesh>
    triSurfaceFrontSphericalTensorField;
typedef Foam::DimensionedField<symmTensor, triSurfaceFrontGeoMesh>
    triSurfaceFrontSymmTensorField;
typedef Foam::DimensionedField<tensor, triSurfaceFrontGeoMesh>
    triSurfaceFrontTensorField;

class triSurfaceFrontPointGeoMesh;

typedef Foam::DimensionedField<label, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointLabelField;
typedef Foam::DimensionedField<scalar, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointScalarField;
typedef Foam::DimensionedField<vector, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointVectorField;
typedef Foam::DimensionedField<sphericalTensor, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointSphericalTensorField;
typedef Foam::DimensionedField<symmTensor, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointSymmTensorField;
typedef Foam::DimensionedField<tensor, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
