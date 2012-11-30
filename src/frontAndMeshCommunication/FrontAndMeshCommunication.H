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
    Foam::FrontAndMeshCommunication

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

Description

SourceFiles
    FrontAndMeshCommunication.C

\*---------------------------------------------------------------------------*/

#ifndef FrontAndMeshCommunication_H
#define FrontAndMeshCommunication_H

#include "tmp.H"
#include "DimensionedField.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "triSurfaceGeoMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace frontTracking
{

/*---------------------------------------------------------------------------*\
                         Class FrontAndMeshCommunication Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh, class Front, class Connectivity>
class FrontAndMeshCommunication
{
    // Private data
    
    const Mesh& mesh_; 
    const Front& front_;

    Connectivity connectivity_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        FrontAndMeshCommunication(const FrontAndMeshCommunication<Mesh, Front, Connectivity>&);

        //- Disallow default bitwise assignment
        void operator=(const FrontAndMeshCommunication<Mesh, Front, Connectivity>&);


public:

    // Constructors

        //- Construct from components
        FrontAndMeshCommunication(const Mesh& mesh, const Front& front);

    // Member Functions 
        void calcCellsToElementsDist(volScalarField& cellsToElementsDist); 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "FrontAndMeshCommunication.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
