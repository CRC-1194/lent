/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "naiveNarrowBandPropagation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

naiveNarrowBandPropagation::naiveNarrowBandPropagation() {}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void naiveNarrowBandPropagation::operator()(volScalarField& Psi)
{
    const fvMesh& mesh = Psi.mesh(); 
    const labelList& own = mesh.owner(); 
    const labelList& nei = mesh.neighbour(); 

    label jumpFace = -1; 

    Psi.time().cpuTimeIncrement(); 
    do 
    {
        jumpFace = -1; 
        forAll (own, faceI)
        {
            if ((Psi[own[faceI]] < 0) && (Psi[nei[faceI]] == GREAT))
            {
                jumpFace = faceI; 
                Psi[nei[faceI]] *= -1; 
            }
            if ((Psi[nei[faceI]] < 0) && (Psi[own[faceI]] == GREAT))
            {
                jumpFace = faceI; 
                Psi[own[faceI]] *= -1;
            }
        }
    } while (jumpFace >= 0);

    Psi.boundaryField().evaluate(); 
    Info << "Enforcing narrow band: " << Psi.time().cpuTimeIncrement() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
