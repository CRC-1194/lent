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
    Foam::FrontAndMeshConnectivity

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
    FrontAndMeshConnectivity.C

\*---------------------------------------------------------------------------*/

#ifndef FrontAndMeshConnectivity_H
#define FrontAndMeshConnectivity_H

#include "typeInfo.H"
#include "MeshObject.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace frontTracking
{

/*---------------------------------------------------------------------------*\
                         Class FrontAndMeshConnectivity Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh, class Front, class ConnectivityData> 
class FrontAndMeshConnectivity
//:
    //public MeshObject<Mesh, 
        //FrontAndMeshConnectivity<Mesh, Front, ConnectivityData> >
{
    // Private data
    
            //- Connectivity maps 
         ConnectivityData elementsToCells_;

public:
        
        //TypeName("frontAndMeshConnectivity");

        typedef typename ConnectivityData::MapType MapType; 

    // Constructors

        //- Construct from components
        FrontAndMeshConnectivity(const Mesh& mesh, const Front& front);

        //- Construct as copy
        //FrontAndMeshConnectivity(const FrontAndMeshConnectivity&);


    //- Destructor
    virtual ~FrontAndMeshConnectivity(); 


    // Member Functions

        // Access
        const typename ConnectivityData::MapType& cellsToElements();

        // Check

        // Edit

        // Write


    // Member Operators

        //void operator=(const FrontAndMeshConnectivity&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FrontAndMeshConnectivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
