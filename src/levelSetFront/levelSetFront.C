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
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "levelSetFront.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frontTracking::levelSetFront::levelSetFront()
:
    triSurface()
{
}

Foam::frontTracking::levelSetFront::levelSetFront(const volScalarField& psi,
                                                  const scalarField& psiPoint)
    : 
        triSurface()
{
    // Reconstruct the interface as a level 0 of the distance field.
    this->reconstruct(psi, psiPoint); 
}


Foam::frontTracking::levelSetFront::levelSetFront(const levelSetFront& rhs)
:
    triSurface(rhs)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::frontTracking::levelSetFront::~levelSetFront()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::frontTracking::levelSetFront::reconstruct
(
    const volScalarField& cellIsoVals, 
    const scalarField& pointIsoVals, 
    const bool regularise, 
    const scalar mergeTol
)
{
    *this = isoSurface(cellIsoVals, pointIsoVals, 0 , regularise, mergeTol);
}

// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

void Foam::frontTracking::levelSetFront::operator=(const isoSurface& rhs)
{
    triSurface::operator=(rhs);
}

// ************************************************************************* //
