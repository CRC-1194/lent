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

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#include "lentMethod.H"

namespace Foam  { 

namespace FrontTracking {

    defineTypeNameAndDebug(lentMethod, 0); 


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType lentMethod::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentMethod::lentMethod(const triSurfaceFront& front, const fvMesh& mesh)
:
    lentControlDict_(
         IOobject(
            "lentSolution", 
            "system", 
            mesh, 
            IOobject::MUST_READ_IF_MODIFIED, 
            IOobject::NO_WRITE
         )
     ),
     distanceFieldCalculatorTmp_(),
     //frontReconstructorTmp_(), 
     //frontVelocityCalculatorTmp_(), 
     //frontMotionSolver_(), 
     heavisideModelTmp_()
{}

lentMethod::lentMethod(const lentMethod& copy)
: 
    lentControlDict_(copy.lentControlDict_),
    distanceFieldCalculatorTmp_(copy.distanceFieldCalculatorTmp_),
    heavisideModelTmp_(copy.heavisideModelTmp_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lentMethod::~lentMethod()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void lentMethod::calcSignedDistanceFields(
    volScalarField& signedDistance, 
    pointScalarField& pointSignedDistance
) const
{
    // Update points-front distance field.
    //distanceFieldCalculatorTmp_->calcPointDistanceField(); 
    // Update cells-front distance field.  
}

void lentMethod::calcHeavisideField(
   volScalarField& heaviside,
   const volScalarField& signedDistance
) const
{
    // Calculate the heaviside field from th 
    // Update cells-front distance field.  
}


void lentMethod::reconstructFront(
    triSurfaceFront& front, 
    const volScalarField& signedDistance, 
    const pointScalarField& pointSignedDistance
) const
{}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void lentMethod::operator=(const lentMethod& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("lentMethod::operator=(const lentMethod&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    distanceFieldCalculatorTmp_ = rhs.distanceFieldCalculatorTmp_; 
    heavisideModelTmp_ = rhs.heavisideModelTmp_; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
