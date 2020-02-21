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

#ifndef triSurfaceFrontFieldsFwd_H
#define triSurfaceFrontFieldsFwd_H

#include "fieldTypes.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

template<class Type, class GeoMesh>
class DimensionedField;

namespace FrontTracking {

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfaceFrontGeoMesh;

typedef DimensionedField<label, triSurfaceFrontGeoMesh>
    triSurfaceFrontLabelField;
typedef DimensionedField<scalar, triSurfaceFrontGeoMesh>
    triSurfaceFrontScalarField;
typedef DimensionedField<vector, triSurfaceFrontGeoMesh>
    triSurfaceFrontVectorField;
typedef DimensionedField<sphericalTensor, triSurfaceFrontGeoMesh>
    triSurfaceFrontSphericalTensorField;
typedef DimensionedField<symmTensor, triSurfaceFrontGeoMesh>
    triSurfaceFrontSymmTensorField;
typedef DimensionedField<tensor, triSurfaceFrontGeoMesh>
    triSurfaceFrontTensorField;

class triSurfaceFrontPointGeoMesh;

typedef DimensionedField<label, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointLabelField;
typedef DimensionedField<scalar, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointScalarField;
typedef DimensionedField<vector, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointVectorField;
typedef DimensionedField<sphericalTensor, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointSphericalTensorField;
typedef DimensionedField<symmTensor, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointSymmTensorField;
typedef DimensionedField<tensor, triSurfaceFrontPointGeoMesh>
    triSurfaceFrontPointTensorField;

class triSurfaceFrontEdgeGeoMesh;

typedef DimensionedField<label, triSurfaceFrontEdgeGeoMesh>
    triSurfaceFrontEdgeLabelField;
typedef DimensionedField<scalar, triSurfaceFrontEdgeGeoMesh>
    triSurfaceFrontEdgeScalarField;
typedef DimensionedField<vector, triSurfaceFrontEdgeGeoMesh>
    triSurfaceFrontEdgeVectorField;
typedef DimensionedField<sphericalTensor, triSurfaceFrontEdgeGeoMesh>
    triSurfaceFrontEdgeSphericalTensorField;
typedef DimensionedField<symmTensor, triSurfaceFrontEdgeGeoMesh>
    triSurfaceFrontEdgeSymmTensorField;
typedef DimensionedField<tensor, triSurfaceFrontEdgeGeoMesh>
    triSurfaceFrontEdgeTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
