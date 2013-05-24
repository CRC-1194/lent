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
    Foam::FvMeshAndFrontConnection

Description

Author
    Tomislav Maric
    tomislav.maric@gmx.com

SourceFiles
    FvMeshAndFrontConnection.C

\*---------------------------------------------------------------------------*/

#ifndef FvMeshAndFrontConnection_H
#define FvMeshAndFrontConnection_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class FvMeshAndFrontConnection Declaration
\*---------------------------------------------------------------------------*/

template<class Front>
class FvMeshAndFrontConnection
{

private: 
    // Private data
    
    const fvMesh& mesh_; 
    const Front& front_;

public:

    // Constructors

        //- Construct from components
        FvMeshAndFrontConnection(const fvMesh& mesh, const Front& front);

        //- Construct from components
        FvMeshAndFrontConnection(const Front& front, const fvMesh& mesh);

    // Member Functions
    
        //- Access
        const fvMesh& mesh() const
        {
            return mesh_; 
        }

        const Front& front() const
        {
            return front_; 
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "FvMeshAndFrontConnection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
