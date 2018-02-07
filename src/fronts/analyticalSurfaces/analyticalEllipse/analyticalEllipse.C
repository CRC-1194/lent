/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::analyticalEllipse

SourceFiles
    analyticalEllipse.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Specialization of the analyticalSurface class for an ellipse.

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#include "analyticalEllipse.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>
#include <assert.h>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalEllipse, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalEllipse, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void analyticalEllipse::ensureValidHalfAxes()
{
    // For this class to work correctly the half axis value of the
    // empty direction can be anything but zero
    forAll(halfAxes_, I)
    {
        if (halfAxes_[I] == 0.0)
        {
            halfAxes_[I] = 1.0;
        }
    }
}

void analyticalEllipse::ensureValidCentre()
{
    // The empty direction component of the centre has to be zero
    // for this class to work properly
    centre_ = projector_&centre_;
}

scalar analyticalEllipse::levelSetValueOf(const point& aPoint) const
{
    // This function resembles
    // f(aPoint) = (x/a)^2 + (y/b)^2 + (z/c)^2 - 1
    scalar levelSetValue = -1.0;

    forAll(aPoint, I)
    {
        levelSetValue += pow(aPoint[I]/halfAxes_[I],2.0);
    }

    return levelSetValue; 
}

vector analyticalEllipse::levelSetGradientAt(const point& aPoint) const
{
    vector levelSetGradient{0,0,0};

    if (mag(aPoint - centre_) < SMALL)
    {
        FatalErrorInFunction
            << "Cannot compute level set gradient for a point which "
            << " coincides with the ellipse centre."
            << abort(FatalError);
    }

    forAll(levelSetGradient, I)
    {
        levelSetGradient[I] = 2.0/pow(halfAxes_[I],2.0)*aPoint[I];
    }

    return levelSetGradient;
}

point analyticalEllipse::moveToReferenceFrame(const point& aPoint) const
{
    return (projector_&aPoint) - centre_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalEllipse::analyticalEllipse(const dictionary& configDict)
:
    analyticalSurface{configDict},
    centre_{configDict.lookup("centre")},
    halfAxes_{configDict.lookup("halfAxes")}
{
    vector emptyDirection = configDict.lookup("emptyDirection");
    projector_ = Identity<scalar>{} - emptyDirection*emptyDirection;
    centre_ = projector_&centre_;
    ensureValidHalfAxes();
    ensureValidCentre();
}

analyticalEllipse::analyticalEllipse(const point& centre, const vector& halfAxes, const vector& emptyDirection)
:
   analyticalSurface{},
   centre_{centre},
   halfAxes_{halfAxes}
{
    projector_ = Identity<scalar>{} - emptyDirection*emptyDirection;
    centre_ = projector_&centre_;
    ensureValidHalfAxes();
    ensureValidCentre();
} 


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar analyticalEllipse::distance(const point& trialPoint) const
{
    return fabs(signedDistance(trialPoint));
}

scalar analyticalEllipse::signedDistance(const point& trialPoint) const
{
    auto copyPoint{trialPoint};
    auto projectedPoint = moveToReferenceFrame(trialPoint);
    auto surfaceProjection = normalProjectionToSurface(copyPoint);

    return sign(levelSetValueOf(projectedPoint))*mag(surfaceProjection - projectedPoint);
}

point analyticalEllipse::normalProjectionToSurface(point& trialPoint) const
{
    auto refPoint = moveToReferenceFrame(trialPoint);

    auto lGrad = levelSetGradientAt(refPoint);

    // Idea to find the normal projection on the ellipse: set up a line
    // using the refPoint and the level set gradient (which is normal to
    // the ellipse) as direction. Intersect this line with the level set
    // description of the ellipse which results in a quadratic equation to
    // be solved
    scalar a = 0.0;
    scalar b = 0.0;
    scalar c = -1.0;

    forAll(refPoint, I)
    {
        a += std::pow(lGrad[I]/halfAxes_[I],2.0);
        b += 2.0*refPoint[I]*lGrad[I]/std::pow(halfAxes_[I],2.0);
        c += std::pow(refPoint[I]/halfAxes_[I],2.0);
    }

    // p-q formula...
    auto p = b/a;
    auto q = c/a;

    // Since the directional vector is normal to the ellipse
    // there must be two real valued solutions for the factor lambda
    // corresponding to the two intersections of the line with the
    // ellipse
    auto lambda1 = -0.5*p + std::sqrt(std::pow(0.5*p,2.0) - q);
    auto lambda2 = -0.5*p - std::sqrt(std::pow(0.5*p,2.0) - q);

    // We want the intersection closer to refPoint, thus the lambda
    // with the smaller absolute value is the one we are looking for
    scalar lambda = 0.0;
    mag(lambda1) < mag(lambda2) ? lambda = lambda1 : lambda = lambda2;

    point ellipsePoint = refPoint + lambda*lGrad;

    assert(mag(levelSetValueOf(ellipsePoint)) < SMALL 
            && "Bug: projected point is not on ellipse");

    // Reverse movement to reference system
    return ellipsePoint + centre_ + 
                ((Identity<scalar>{} - projector_)&trialPoint);
}

vector analyticalEllipse::normalToPoint(const point& trialPoint) const
{
    auto refPoint = moveToReferenceFrame(trialPoint);
    auto gradient = levelSetGradientAt(refPoint);
    auto levelSetValue = levelSetValueOf(refPoint);

    return sign(levelSetValue)*gradient/mag(gradient);
}

point analyticalEllipse::intersection(const point& pointA,
                                     const point& pointB) const
{
    notImplemented("analyticalEllipsoid::intersection(...)");
    return point{0,0,0};
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
analyticalEllipse& analyticalEllipse::operator=(const analyticalEllipse& rhs)
{
    if (this != &rhs)
    {
        centre_ = rhs.centre_;
        halfAxes_ = rhs.halfAxes_;
        projector_ = rhs.projector_;
    }

    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
