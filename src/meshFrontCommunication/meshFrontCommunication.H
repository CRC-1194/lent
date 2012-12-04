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
    Foam::meshFrontCommunication

Description

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany


SourceFiles
    meshFrontCommunication.C

\*---------------------------------------------------------------------------*/

#ifndef meshFrontCommunication_H
#define meshFrontCommunication_H

#include "levelSetFront.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace frontTracking
{

/*---------------------------------------------------------------------------*\
                         Class meshFrontCommunication Declaration
\*---------------------------------------------------------------------------*/

class meshFrontCommunication
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        //meshFrontCommunication(const meshFrontCommunication&);

        //- Disallow default bitwise assignment
        //void operator=(const meshFrontCommunication&);


public:

    // Static data members

    // Constructors

        //- Construct null
        //meshFrontCommunication();

        //- Construct from components
        meshFrontCommunication(const fvMesh& mesh, const levelSetFront& front);

        //- Construct from Istream
        //meshFrontCommunication(Istream&);

        //- Construct as copy
        //meshFrontCommunication(const meshFrontCommunication&);


    //- Destructor
    //~meshFrontCommunication();


    // Member Functions

        // Calculate 
        void calcDistanceField (volScalarField& Psi); 

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        //void operator=(const meshFrontCommunication&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
