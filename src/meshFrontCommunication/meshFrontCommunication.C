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

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt

\*---------------------------------------------------------------------------*/

#include "meshFrontCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::frontTracking::meshFrontCommunication::meshFrontCommunication (
    const fvMesh& mesh, 
    const levelSetFront& front
)
{
    Pout << " Foam::frontTracking::meshFrontCommunication::meshFrontCommunication ( "
        << "\n    const fvMesh& mesh, " 
        << "\n    const levelSetFront& front " << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//Foam::frontTracking::meshFrontCommunication::~meshFrontCommunication()
//{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::frontTracking::meshFrontCommunication::calcDistanceField(volScalarField& psi)
{
    Pout << "Foam::frontTracking::meshFrontCommunication::calcDistanceField(volScalarField& psi)" 
        << endl;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//void Foam::frontTracking::meshFrontCommunication::operator=(const meshFrontCommunication& rhs)
//{
    //// Check for assignment to self
    //if (this == &rhs)
    //{
        //FatalErrorIn("Foam::frontTracking::meshFrontCommunication::operator=(const Foam::frontTracking::meshFrontCommunication&)")
            //<< "Attempted assignment to self"
            //<< abort(FatalError);
    //}
//}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
