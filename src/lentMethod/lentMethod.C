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
#include "error.H"

namespace Foam  { 

namespace FrontTracking {

    defineTypeNameAndDebug(lentMethod, 0); 


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType lentMethod::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentMethod::lentMethod(
    const triSurfaceFront& front, 
    const fvMesh& mesh,
    word dictName 
)
:
    regIOobject(
        IOobject(
            "lentMethod", 
            mesh.time().timeName(), 
            mesh.time(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        )
    ),
    frontTmp_(front), 
    lentControlDict_(
         IOobject(
            dictName, 
            "system", 
            mesh.time(), 
            IOobject::MUST_READ_IF_MODIFIED, 
            IOobject::AUTO_WRITE
         )
     ),
     lentDistanceFieldCalculatorTmp_(
        lentDistanceFieldCalculator::New(
            lentControlDict_.subDict("distanceCalculator")
        )
     ),
     //frontReconstructorTmp_(), 
     //frontVelocityCalculatorTmp_(), 
     //frontMotionSolver_(), 
     heavisideModelTmp_(
         heavisideModel::New(
            lentControlDict_.lookup("heavisideModel"), 
            lentControlDict_.subDict("heavisideModel")
         )
     ) 
{
}

lentMethod::lentMethod(const lentMethod& copy)
: 
    regIOobject(copy),
    lentControlDict_(copy.lentControlDict_),
    //lentDistanceFieldCalculatorTmp_(copy.lentDistanceFieldCalculatorTmp_),
    heavisideModelTmp_(copy.heavisideModelTmp_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lentMethod::~lentMethod()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void lentMethod::calcSearchDistances(
    volScalarField& searchDistanceSqr, 
    pointScalarField& pointSearchDistanceSqr
) 
{
    lentDistanceFieldCalculator& distanceCalc = lentDistanceFieldCalculatorTmp_(); 

    distanceCalc.calcCellSearchDistance(searchDistanceSqr); 
    distanceCalc.calcPointSearchDistance(pointSearchDistanceSqr, searchDistanceSqr); 
}

void lentMethod::calcSignedDistances(
    volScalarField& signedDistance, 
    pointScalarField& pointSignedDistance
) 
{
    // Update points-front distance field.
    //lentDistanceFieldCalculatorTmp_->calcPointsToFrontDistanceField(); 
    // Update cells-front distance field.  
    //lentDistanceFieldCalculatorTmp_->cellsToFrontDistanceField(); 
}

void lentMethod::calcHeaviside(
   volScalarField& heaviside,
   const volScalarField& signedDistance
) const
{
    heavisideModelTmp_->calcHeaviside(
        heaviside, 
        signedDistance
    ); 
}

void lentMethod::reconstructFront(
    triSurfaceFront& front, 
    const volScalarField& signedDistance
) const
{}

//void lentMethod::setHeavisideModel(word name) 
//{
    //heavisideModelTmp_ = heavisideModel::New(name); 
//}

bool lentMethod::writeData(Ostream& os) const
{
    FatalErrorIn("lentMethod::writeData(Ostream& os)")
    << "lentMethod is not supposed to be written "
    << "regIOobject inherited to allow registry queries." << endl;

    return false; 
}

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

    lentDistanceFieldCalculatorTmp_ = rhs.lentDistanceFieldCalculatorTmp_; 
    heavisideModelTmp_ = rhs.heavisideModelTmp_; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //