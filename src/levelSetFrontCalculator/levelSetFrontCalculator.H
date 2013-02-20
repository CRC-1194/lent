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
    Foam::levelSetFrontCalculator

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
    levelSetFrontCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef levelSetFrontCalculator_H
#define levelSetFrontCalculator_H

#include "levelSetFront.H"
#include "fvMesh.H"
#include "levelSetFrontAndMeshConnection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace levelSetFrontTracking 
{

/*---------------------------------------------------------------------------*\
                         Class levelSetFrontCalculator Declaration
\*---------------------------------------------------------------------------*/

class levelSetFrontCalculator
{
    // Private data
    levelSetFrontAndMeshConnection connectivity_; 

    // Private Member Functions

    void resetDistanceField(volScalarField& psi); 
    scalar pointToElementDistance(const point& inputPoint, label elementI); 

        //- Disallow default bitwise copy construct
        //levelSetFrontCalculator(const levelSetFrontCalculator&);

        //- Disallow default bitwise assignment
        //void operator=(const levelSetFrontCalculator&);


public:

    // Static data members

    // Constructors

        //- Construct null
        //levelSetFrontCalculator();

        //- Construct from components
        levelSetFrontCalculator(const fvMesh& mesh, const levelSetFront& front);

        //- Construct from Istream
        //levelSetFrontCalculator(Istream&);

        //- Construct as copy
        //levelSetFrontCalculator(const levelSetFrontCalculator&);


    //- Destructor
    //~levelSetFrontCalculator();


    // Member Functions

        // Calculate 
        void calcDistanceField (volScalarField& Psi); 

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        //void operator=(const levelSetFrontCalculator&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace levelSetFrontTracking 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
