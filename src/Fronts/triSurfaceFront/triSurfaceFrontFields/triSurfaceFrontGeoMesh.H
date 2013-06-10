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
    The triSurfaceFront GeoMesh (for holding fields).

    Similar to the volGeoMesh used for the Finite Volume discretization.

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceFrontGeoMesh_H
#define triSurfaceFrontGeoMesh_H

#include "GeoMesh.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfaceFrontGeoMesh
:
    public GeoMesh<triSurfaceFront>
{

public:

    // Constructors

        //- Construct from triSurfaceFront reference
        explicit triSurfaceFrontGeoMesh(const triSurfaceFront& mesh)
        :
            GeoMesh<triSurfaceFront>(mesh)
        {}


    // Member Functions

        //- Return size
        static label size(const triSurfaceFront& mesh)
        {
            return mesh.size();
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
