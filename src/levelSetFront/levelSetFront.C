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
    Germany

\*---------------------------------------------------------------------------*/

#include "levelSetFront.H"
#include "volPointInterpolation.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * //

namespace Foam
{
    namespace frontTracking
    {
        defineTypeNameAndDebug(levelSetFront, 0);
    }
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void Foam::frontTracking::levelSetFront::computeIsoSurface (
    const volScalarField& cellsToElementsDist, 
    const scalarField& pointsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    *this = isoSurface (
        cellsToElementsDist, 
        pointsToElementsDist, 
        0, 
        regularise,
        mergeTol
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frontTracking::levelSetFront::levelSetFront(const IOobject& io)
:
    triSurfaceMesh(io), 
    moving_(false), 
    changing_(false)
{
}

//Foam::frontTracking::levelSetFront::levelSetFront(const volScalarField& psi,
                                                  //const scalarField& psiPoint)
    //: 
        //triSurfaceMesh(), 
        //moving_(false), 
        //changing_(false)
//{
    //// Reconstruct the interface as a level 0 of the distance field.
    //computeIsoSurface(psi, psiPoint); 
//}


//Foam::frontTracking::levelSetFront::levelSetFront(const levelSetFront& rhs)
//:
    //triSurfaceMesh(rhs),
    //moving_(rhs.moving_),
    //changing_(rhs.changing_)
//{

//}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//Foam::frontTracking::levelSetFront::~levelSetFront()
//{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::frontTracking::levelSetFront::moving() const
{
    return moving_;
}
void Foam::frontTracking::levelSetFront::moving(bool b)
{
    moving_ = b;
}

bool Foam::frontTracking::levelSetFront::changing() const
{
    return changing_;
}
void Foam::frontTracking::levelSetFront::changing(bool b)
{
    changing_ = b;
}

void Foam::frontTracking::levelSetFront::reconstruct (
    const volScalarField& cellsToElementsDist, 
    const scalarField& pointsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    computeIsoSurface (
        cellsToElementsDist, 
        pointsToElementsDist, 
        regularise, mergeTol
    );

    changing(true);
}

void Foam::frontTracking::levelSetFront::reconstruct (
    const volScalarField& cellsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    const fvMesh& mesh = cellsToElementsDist.mesh(); 
    volPointInterpolation pInter (mesh);

    tmp<pointScalarField> pointsToElementsDistTmp = pInter.interpolate (
        cellsToElementsDist
    );

    const pointScalarField& pointsToElementsDist = pointsToElementsDistTmp();

    computeIsoSurface (
        cellsToElementsDist, 
        pointsToElementsDist, 
        regularise, mergeTol
    );

    changing(true);
}

// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

void Foam::frontTracking::levelSetFront::operator=(const isoSurface& rhs)
{
    triSurface::operator=(rhs);
}

// ************************************************************************* //
