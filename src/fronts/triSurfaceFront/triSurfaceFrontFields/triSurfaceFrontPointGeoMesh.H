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

Class
    Foam::triSurfaceFrontGeoMesh

Description
    The triSurfaceFront  point GeoMesh (for holding vertex fields).

    Similar to the volMesh used for the Finite Volume discretization.

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceFrontPointGeoMesh_H
#define triSurfaceFrontPointGeoMesh_H

#include "GeoMesh.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfaceFrontPointGeoMesh
:
    public GeoMesh<triSurfaceFront>
{

public:

    // Constructors

        //- Construct from triSurfaceFront reference
        explicit triSurfaceFrontPointGeoMesh(const triSurfaceFront& mesh)
        :
            GeoMesh<triSurfaceFront>(mesh)
        {}


    // Member Functions

        //- Return size
        static label size(const triSurfaceFront& mesh)
        {
            const triSurface& surface = mesh;
            return surface.points().size();
        }

        //- Return size
        label size() const
        {
            return size(mesh_);
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
