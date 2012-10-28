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
    Foam::frontAndMeshConnectivity

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

Description

SourceFiles
    frontAndMeshConnectivity.C

\*---------------------------------------------------------------------------*/

#ifndef frontAndMeshConnectivity_H
#define frontAndMeshConnectivity_H

#include "typeInfo.H"
#include "MeshObject.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace frontTracking
{

/*---------------------------------------------------------------------------*\
                         Class frontAndMeshConnectivity Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh> 
class frontAndMeshConnectivity
:
    public MeshObject<Mesh, frontAndMeshConnectivity<Mesh> >
{

    // Private data
        //- Front vertex to cell  
        //- Front element to cell

        //- Cell to front vertices

        //- Cell to front elements

    // Private Member Functions


    // Dissallow
        //- Empty construction
        //frontAndMeshConnectivity<Mesh>();

public:
        

    // Constructors

        TypeName("frontAndMeshConnectivity");


        //- Construct from components
        frontAndMeshConnectivity(const Mesh& mesh, 
                                 const levelSetFront& front);

        //- Construct as copy
        //frontAndMeshConnectivity(const frontAndMeshConnectivity&);


    // Selectors

        //- Select null constructed
        //static autoPtr<frontAndMeshConnectivity> New();

    //- Destructor
    ~frontAndMeshConnectivity();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        //void operator=(const frontAndMeshConnectivity&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "frontAndMeshConnectivityI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
